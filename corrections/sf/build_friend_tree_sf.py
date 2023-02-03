from sys import base_prefix
import ROOT
import argparse
import yaml
import os
import time
import glob
import shutil
from tqdm import tqdm
from multiprocessing import Pool, current_process, RLock
from collections import OrderedDict

def base_filename(path):
    return path.split("/")[-3]

def job_wrapper(args):
    # print(args)
    return friend_producer(*args)

def friend_producer(rfile, scalefactors):
    output_path = rfile.replace("ntuples_xsec", "ntuples_xsec_sf_EraC")
    # output_path = rfile.replace("ntuples", "friends/sf")

    if os.path.exists(output_path):
        print(f"friend_producer: {output_path} exists -> skip")
        return

    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=False)

    rdf = ROOT.RDataFrame("ntuple", rfile)
    original_cols = [str(col) for col in rdf.GetColumnNames()]

    isdilepton = ("/mm/" in rfile or "/ee/" in rfile)

    out_columns = []

    for type, info in scalefactors.items():
        for key in info.keys():
            if key == "name":
                continue
            if type == "trg":
                val_name = f"val_{type}_{key}_1"
                val_str = f"(float){info[key]}->GetBinContent(\
                    {info[key]}->GetXaxis()->FindBin(q_1),\
                    {info[key]}->GetYaxis()->FindBin(TMath::Max((float)(-2.4+ 1.e-3), TMath::Min((float)(2.4 - 1.e-3), (float)(eta_1)))),\
                    {info[key]}->GetZaxis()->FindBin(TMath::Max((float)(25. + 1.e-3), TMath::Min((float)(120. - 1.e-3), (float)pt_1)))\
                )"
                rdf = rdf.Define(val_name, val_str)
                out_columns.append(val_name)

                err_name = f"err_{type}_{key}_1"
                err_str = f"(float){info[key]}->GetBinError(\
                    {info[key]}->GetXaxis()->FindBin(q_1),\
                    {info[key]}->GetYaxis()->FindBin(TMath::Max((float)(-2.4+ 1.e-3), TMath::Min((float)(2.4 - 1.e-3), (float)(eta_1)))),\
                    {info[key]}->GetZaxis()->FindBin(TMath::Max((float)(25. + 1.e-3), TMath::Min((float)(120. - 1.e-3), (float)pt_1)))\
                )"
                rdf = rdf.Define(err_name, err_str)
                out_columns.append(err_name)
            else:
                val_name = f"val_{type}_{key}_1"
                val_str = f"(float){info[key]}->GetBinContent(\
                    {info[key]}->GetXaxis()->FindBin(TMath::Max((float)(-2.4+ 1.e-3), TMath::Min((float)(2.4 - 1.e-3), (float)(eta_1)))),\
                    {info[key]}->GetYaxis()->FindBin(TMath::Max((float)(25. + 1.e-3), TMath::Min((float)(120. - 1.e-3), (float)pt_1)))\
                )"
                rdf = rdf.Define(val_name, val_str)
                out_columns.append(val_name)

                err_name = f"err_{type}_{key}_1"
                err_str = f"(float){info[key]}->GetBinError(\
                    {info[key]}->GetXaxis()->FindBin(TMath::Max((float)(-2.4+ 1.e-3), TMath::Min((float)(2.4 - 1.e-3), (float)(eta_1)))),\
                    {info[key]}->GetYaxis()->FindBin(TMath::Max((float)(25. + 1.e-3), TMath::Min((float)(120. - 1.e-3), (float)pt_1)))\
                )"
                rdf = rdf.Define(err_name, err_str)
                out_columns.append(err_name)

            rdf = rdf.Define(val_name+"_up", f"{val_name} + {err_name}")
            rdf = rdf.Define(val_name+"_dn", f"{val_name} - {err_name}")
            out_columns.append(val_name+"_up")
            out_columns.append(val_name+"_dn")

            if isdilepton:
                if type == "trg":
                    val_name = f"val_{type}_{key}_2"
                    val_str = f"(float){info[key]}->GetBinContent(\
                        {info[key]}->GetXaxis()->FindBin(q_2),\
                        {info[key]}->GetYaxis()->FindBin(TMath::Max((float)(-2.4+ 1.e-3), TMath::Min((float)(2.4 - 1.e-3), (float)(eta_2)))),\
                        {info[key]}->GetZaxis()->FindBin(TMath::Max((float)(25. + 1.e-3), TMath::Min((float)(120. - 1.e-3), (float)pt_2)))\
                    )"
                    rdf = rdf.Define(val_name, val_str)
                    out_columns.append(val_name)

                    err_name = f"err_{type}_{key}_2"
                    err_str = f"(float){info[key]}->GetBinError(\
                        {info[key]}->GetXaxis()->FindBin(q_2),\
                        {info[key]}->GetYaxis()->FindBin(TMath::Max((float)(-2.4+ 1.e-3), TMath::Min((float)(2.4 - 1.e-3), (float)(eta_2)))),\
                        {info[key]}->GetZaxis()->FindBin(TMath::Max((float)(25. + 1.e-3), TMath::Min((float)(120. - 1.e-3), (float)pt_2)))\
                    )"
                    rdf = rdf.Define(err_name, err_str)
                    out_columns.append(err_name)
                else:
                    val_name = f"val_{type}_{key}_2"
                    val_str = f"(float){info[key]}->GetBinContent(\
                        {info[key]}->GetXaxis()->FindBin(TMath::Max((float)(-2.4+ 1.e-3), TMath::Min((float)(2.4 - 1.e-3), (float)(eta_2)))),\
                        {info[key]}->GetYaxis()->FindBin(TMath::Max((float)(25. + 1.e-3), TMath::Min((float)(120. - 1.e-3), (float)pt_2)))\
                    )"
                    rdf = rdf.Define(val_name, val_str)
                    out_columns.append(val_name)

                    err_name = f"err_{type}_{key}_2"
                    err_str = f"(float){info[key]}->GetBinError(\
                        {info[key]}->GetXaxis()->FindBin(TMath::Max((float)(-2.4+ 1.e-3), TMath::Min((float)(2.4 - 1.e-3), (float)(eta_2)))),\
                        {info[key]}->GetYaxis()->FindBin(TMath::Max((float)(25. + 1.e-3), TMath::Min((float)(120. - 1.e-3), (float)pt_2)))\
                    )"
                    rdf = rdf.Define(err_name, err_str)
                    out_columns.append(err_name)

                rdf = rdf.Define(val_name+"_up", f"{val_name} + {err_name}")
                rdf = rdf.Define(val_name+"_dn", f"{val_name} - {err_name}")
                out_columns.append(val_name+"_up")
                out_columns.append(val_name+"_dn")

        sf_name = f"sf_{type}"
        if isdilepton:
            if type == "trg":
                eff_eff_dt    = f"(1. - (1. - (val_{type}_mc_1 * val_{type}_sf_1))   *(1. - (val_{type}_mc_2 * val_{type}_sf_2)))"
                eff_eff_dt_up = f"(1. - (1. - (val_{type}_mc_1 * val_{type}_sf_1_up))*(1. - (val_{type}_mc_2 * val_{type}_sf_2_up)))"
                eff_eff_dt_dn = f"(1. - (1. - (val_{type}_mc_1 * val_{type}_sf_1_dn))*(1. - (val_{type}_mc_2 * val_{type}_sf_2_dn)))"
                eff_eff_mc = f"(1. - (1. - val_{type}_mc_1)*(1. - val_{type}_mc_2))"
                rdf = rdf.Define(sf_name, f"(float)({eff_eff_dt}/{eff_eff_mc})")
                rdf = rdf.Define(sf_name+"_up", f"(float)({eff_eff_dt_up}/{eff_eff_mc})")
                rdf = rdf.Define(sf_name+"_dn", f"(float)({eff_eff_dt_dn}/{eff_eff_mc})")
            else:
                rdf = rdf.Define(sf_name, f"(float)(val_{type}_sf_1)*(val_{type}_sf_2)")
                rdf = rdf.Define(sf_name+"_up", f"(float)(val_{type}_sf_1_up)*(val_{type}_sf_2_up)")
                rdf = rdf.Define(sf_name+"_dn", f"(float)(val_{type}_sf_1_dn)*(val_{type}_sf_2_dn)")
        else:
            rdf = rdf.Define(sf_name, f"(float)(val_{type}_sf_1)")
            rdf = rdf.Define(sf_name+"_up", f"(float)(val_{type}_sf_1_up)")
            rdf = rdf.Define(sf_name+"_dn", f"(float)(val_{type}_sf_1_dn)")
        out_columns.append(sf_name)
        out_columns.append(sf_name+"_up")
        out_columns.append(sf_name+"_dn")


    rdf.Snapshot(
        "ntuple",
        output_path,
        original_cols + out_columns
    )

def generate_friend_trees(scalefactors, ntuples, nthreads):
    arguments = [(ntuple, scalefactors) for ntuple in ntuples]

    pool = Pool(nthreads, initargs=(RLock(),), initializer=tqdm.set_lock)
    for _ in tqdm(
        pool.imap_unordered(job_wrapper, arguments),
        total=len(arguments),
        desc="Total progess",
        # position=nthreads + 1,
        dynamic_ncols=True,
        leave=True,
    ):
        pass


if __name__ == "__main__":
    ROOT.TH1.AddDirectory(ROOT.kFALSE)
    ROOT.TH2.AddDirectory(ROOT.kFALSE)
    ROOT.TH3.AddDirectory(ROOT.kFALSE)
    ROOT.gROOT.SetBatch(True)
    ROOT.ROOT.EnableImplicitMT(24)

    # base_path = "/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec/2022/*/mm*/*.root"
    base_path = "/work/jdriesch/earlyrun3/samples/Run3V04/ntuples_xsec/2022/*/mm*/*.root"

    scalefactors = OrderedDict({
        "trk": {
            "name": "NUM_GlobalMuons_DEN_StandAlone_EraC_etasingle_pt",
        },
        "sta": {
            "name": "NUM_StandAlone_DEN_genTracks_EraC_eta_pt",
        },
        "glb": {
            "name": "NUM_GlobalMuons_DEN_genTracks_EraC_eta_pt",
        },
        "id": {
            "name": "NUM_TightID_DEN_GlobalMuons_EraC_eta_pt",
        },
        "iso": {
            "name": "NUM_PFIsoTight_DEN_TightID_EraC_eta_pt",
        },
        "sel": {
            "name": "NUM_TightID_and_PFIsoTight_DEN_genTracks_EraC_eta_pt",
        },
        "trg": {
            "name": "NUM_IsoMu24_DEN_TightID_and_PFIsoTight_EraC_charge_eta_pt",
        },
    })


    ROOT.gROOT.ProcessLine('TFile* f_eff = TFile::Open("data/Efficiencies_muon_Z_Run2022_EraC.root");')
    for type, info in scalefactors.items():
        h_dt_name = f"h_eff_{type}_dt"
        h_mc_name = f"h_eff_{type}_mc"
        h_sf_name = f"h_eff_{type}_sf"
        if type == "trg":
            ROOT.gROOT.ProcessLine(f'TH3F* {h_dt_name} = (TH3F*)f_eff->Get("{info["name"]}_efficiencyData");')
            ROOT.gROOT.ProcessLine(f'TH3F* {h_mc_name} = (TH3F*)f_eff->Get("{info["name"]}_efficiencyMC");')
            ROOT.gROOT.ProcessLine(f'TH3F* {h_sf_name} = (TH3F*)f_eff->Get("{info["name"]}_syst");')
        else:
            ROOT.gROOT.ProcessLine(f'TH2F* {h_dt_name} = (TH2F*)f_eff->Get("{info["name"]}_efficiencyData");')
            ROOT.gROOT.ProcessLine(f'TH2F* {h_mc_name} = (TH2F*)f_eff->Get("{info["name"]}_efficiencyMC");')
            ROOT.gROOT.ProcessLine(f'TH2F* {h_sf_name} = (TH2F*)f_eff->Get("{info["name"]}_syst");')
        scalefactors[type]['dt'] = h_dt_name
        scalefactors[type]['mc'] = h_mc_name
        scalefactors[type]['sf'] = h_sf_name
    ROOT.gROOT.ProcessLine('f_eff->Close();')

    ntuples = glob.glob(base_path)
    ntuples_wo_data = ntuples.copy()
    for ntuple in ntuples:
        if ("Run2022D" in ntuple) or ("EGamma" in ntuple):
            ntuples_wo_data.remove(str(ntuple))
    nthreads = 64
    if nthreads > len(ntuples_wo_data):
        nthreads = len(ntuples_wo_data)
    print("# files: ", len(ntuples_wo_data))
    generate_friend_trees(scalefactors, ntuples_wo_data, nthreads)
