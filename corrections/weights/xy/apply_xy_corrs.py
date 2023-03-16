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

def apply_corrections(f, corr_file):
    output_path = f.replace("ntuples_xsec_sf_scaleres_ptuncorr_pu_EraC", "friend_xsec_sf_scaleres_ptuncorr_EraC_xycorr")    
    if os.path.isfile(output_path):
        f_tmp = None
        try:
            f_tmp = ROOT.TFile(output_path)
        except OSError:
            os.remove(output_path)
        if f_tmp and not f_tmp.IsZombie():
            return 1  
    
    with open(corr_file) as cf:
        corr_dict = json.load(cf)

    chain = ROOT.TChain("ntuple")
    chain.Add(f)
    friend = ROOT.TChain("ntuple")
    friend.Add(f.replace("ntuples_xsec_sf_scaleres_ptuncorr_pu_EraC", "friend_xsec_sf_scaleres_ptuncorr_lepunc_pu_EraC_met_corr"))

    chain.AddFriend(friend)

    rdf = ROOT.RDataFrame(chain)
    
    # check if data
    is_data = (rdf.Sum("is_data").GetValue()>0)
    # is_sigw = ("WtoLNu" in f and "mmet" in f)

    rdf = rdf.Define("pfmet_corr_x", "pfmet_corr * cos(pfmetphi_corr)")
    rdf = rdf.Define("pfmet_corr_y", "pfmet_corr * sin(pfmetphi_corr)")

    if is_data: 
        ch = 'data'
    else:
        ch = 'mc'

    rdf = rdf.Define("pfmet_xycorr_x", f"pfmet_corr_x - (({corr_dict[ch+'_x']['m']}) * npvGood + ({corr_dict[ch+'_x']['c']}))")
    rdf = rdf.Define("pfmet_xycorr_y", f"pfmet_corr_y - (({corr_dict[ch+'_y']['m']}) * npvGood + ({corr_dict[ch+'_y']['c']}))")

    rdf = rdf.Define("pfmet_xycorr", "sqrt(pfmet_xycorr_x * pfmet_xycorr_x + pfmet_xycorr_y * pfmet_xycorr_y)")
    rdf = rdf.Define("pfmetphi_xycorr", "atan2(pfmet_xycorr_y, pfmet_xycorr_x)")
        
    outdir = output_path.replace(output_path.split('/')[-1], "")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    quants = ["pfmet_xycorr", "pfmetphi_xycorr"]

    rdf.Snapshot("ntuple", output_path, quants)

#"""
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
#"""

if __name__=='__main__':
    # ROOT.gROOT.SetBatch(ROOT.kTRUE)
    # ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    # ROOT.ROOT.EnableImplicitMT(16)
    ROOT.gROOT.SetBatch(True)


    # Load ntuples
    base_path = "/storage/9/jdriesch/earlyrun3/samples/Run3V06/ntuples_xsec_sf_scaleres_ptuncorr_pu_EraC/2022/*/*/*.root"
    # base_path = "/storage/9/jdriesch/earlyrun3/samples/Run3V06/friend_xsec_sf_scaleres_ptuncorr_lepunc_pu_EraC_met_corr/2022/*/*/*.root"
    ntuples = glob.glob(base_path)

    # Load weight files
    corr_file = 'corr.json'

    hrange = {
        'pfmetphi_corr': [-np.pi, np.pi, 100],
    }

    nthreads = 16
    arguments = [(ntuple, corr_file) for ntuple in ntuples]
    # print(arguments[0])
    # apply_corrections(ntuples[0], puweights, hrange)

    generate_files(arguments, nthreads)
    #for n in tqdm(ntuples):
    #    apply_corrections(n, x, mz_mc, mz_dt, pt_sf, mz_res_mc, mz_res_dt)
