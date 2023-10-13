from ast import arg
import ROOT
import os
from tqdm import tqdm
import numpy as np
import array as a
import json
import glob
from multiprocessing import Pool, RLock


def job_wrapper(args):
    return apply_corrections(*args)

def apply_corrections(f, wfile, hrange):
    output_path = f.replace("ntuples", "friends/pu")
    
    if os.path.isfile(output_path):
        f_tmp = None
        try:
            f_tmp = ROOT.TFile(output_path)
        except OSError:
            os.remove(output_path)
        if f_tmp and not f_tmp.IsZombie():
            return 1  
    
    
    wf = ROOT.TFile(wfile, "READ")
    wdict = {}
    for var in hrange.keys():
        wdict[str(var)] = wf.Get(var+"weight")
        wdict[var].SetDirectory(ROOT.nullptr)
    wf.Close()


    rdf = ROOT.RDataFrame('ntuple', f)

    # check if data
    is_data = (rdf.Sum("is_data").GetValue()>0) 
    
    rdf = rdf.Define('puweight', '1.')
    rdf = rdf.Define('puweightUp', '1.')
    rdf = rdf.Define('puweightDn', '1.')

    outdir = output_path.replace(output_path.split('/')[-1], "")
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    rdf.Snapshot("ntuple", output_path, ['puweight', 'puweightUp', 'puweightDn'])


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
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    #ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.gROOT.SetBatch(True)


    # Load ntuples
    base_path = "/ceph/jdriesch/CROWN_samples/Run3V07/ntuples/2022/EWK*/mm*/*.root"
    ntuples = glob.glob(base_path)

    # Load weight files
    rfile_w = 'weights.root'

    hrange = {
        'npvGood': [0, 60, 60],
        'rhoFastjetCentralChargedPileUp': [0, 50, 50],
        'rhoFastjetCentralCalo': [0, 35, 35]
    }    
    nthreads = 16
    arguments = [(ntuple, rfile_w, hrange) for ntuple in ntuples]
    for a in arguments:
        job_wrapper(a)
