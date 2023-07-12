from ast import arg
import ROOT
import os
from tqdm import tqdm
import numpy as np
import array as a
import yaml
import glob
from multiprocessing import Pool, RLock


def job_wrapper(args):
    return apply_corrections(*args)

def apply_corrections(f, corr_dict):
    output_path = f.replace("ntuples", "friends/xy")
    if os.path.isfile(output_path):
        f_tmp = None
        try:
            f_tmp = ROOT.TFile(output_path)
        except OSError:
            os.remove(output_path)
        if f_tmp and not f_tmp.IsZombie():
            return 1  
    
    chain = ROOT.TChain("ntuple")
    chain.Add(f)

    rdf = ROOT.RDataFrame(chain)
    
    # check if data
    is_data = (rdf.Sum("is_data").GetValue()>0)

    rdf = rdf.Define("pfmet_x", "pfmet_uncorrected * cos(pfmetphi_uncorrected)")
    rdf = rdf.Define("pfmet_y", "pfmet_uncorrected * sin(pfmetphi_uncorrected)")

    if is_data: 
        ch = 'data'
    else:
        ch = 'mc'
    
    with open(corr_dict[ch], 'r') as cf:
        corr_yaml = yaml.load(cf, Loader=yaml.Loader)

    rdf = rdf.Define("pfmet_xycorr_x", f"pfmet_x - (({corr_yaml['_x']['m']}) * npvGood + ({corr_yaml['_x']['c']}))")
    rdf = rdf.Define("pfmet_xycorr_y", f"pfmet_y - (({corr_yaml['_y']['m']}) * npvGood + ({corr_yaml['_y']['c']}))")

    rdf = rdf.Define("pfmet_xycorr", "sqrt(pfmet_xycorr_x * pfmet_xycorr_x + pfmet_xycorr_y * pfmet_xycorr_y)")
    rdf = rdf.Define("pfmetphi_xycorr", "atan2(pfmet_xycorr_y, pfmet_xycorr_x)")
        
    outdir = output_path.replace(output_path.split('/')[-1], "")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    quants = ["pfmet_xycorr", "pfmetphi_xycorr"]

    rdf.Snapshot("ntuple", output_path, quants)


def generate_files(arguments, nthreads):
    pool = Pool(nthreads, initargs=(RLock(),), initializer=tqdm.set_lock)
    for _ in tqdm(
        pool.imap_unordered(job_wrapper, arguments),
        total=len(arguments),
        desc="Total progess",
        dynamic_ncols=True,
        leave=True,
        ):
        pass


if __name__=='__main__':
    ROOT.gROOT.SetBatch(True)


    # Load ntuples
    base_path = "/ceph/jdriesch/CROWN_samples/Run3V07/ntuples/2022/*/mm*/*.root"
    ntuples = glob.glob(base_path)

    # Load weight files
    corr_dict = {
        'data': 'corrections/2022C/MET.yaml',
        'mc': 'corrections/Run3Summer22NanoAODv11-126X/MET.yaml',
    }

    nthreads = 16
    arguments = [(ntuple, corr_dict) for ntuple in ntuples]

    generate_files(arguments, nthreads)
