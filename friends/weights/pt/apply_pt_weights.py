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
    output_path = f.replace("ntuples_xsec_sf_scaleres_ptuncorr_pu_EraC", "friend_xsec_sf_scaleres_ptuncorr_EraC_ptweights")
    
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
    # is_data = (rdf.Sum("is_data").GetValue()>0) 
    is_sigw = ("WtoLNu" in f and "mmet" in f)

    quants = []
        
    for var in hrange.keys():
        w = var+'weight'
        bins = np.linspace(hrange[var][0], hrange[var][1], hrange[var][2]+1)

        rdf = rdf.Define(w,'1')

        if is_sigw:
            rdf = rdf.Redefine(w,str(wdict[var].GetBinContent(hrange[var][2]+1))) # also correct overflow bin

            for i in range(hrange[var][2]):
                rdf = rdf.Redefine(
                    w,
                    'double w;'\
                    f'if ({var} >= {bins[i]} && {var} < {bins[i+1]})'\
                    f'w = {wdict[var].GetBinContent(i+1)};'\
                    f'else w = {w};'\
                    'return w;'
                )
        rdf = rdf.Define(f"{w}Up", f"{w} + .3*({w}-1)")
        rdf = rdf.Define(f"{w}Dn", f"{w} - .3*({w}-1)")

        for updn in ["", "Up", "Dn"]:
            mean = rdf.Mean(w+updn).GetValue()
            rdf = rdf.Redefine(w+updn, f"{w}{updn}/{mean}")
            quants.append(w+updn)
        
    outdir = output_path.replace(output_path.split('/')[-1], "")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

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
    ntuples = glob.glob(base_path)

    # Load weight files
    rfile_w = 'weights.root'

    hrange = {
        'bosonpt': [0, 100, 100],
    }

    nthreads = 16
    arguments = [(ntuple, rfile_w, hrange) for ntuple in ntuples]
    # print(arguments[0])
    # apply_corrections(ntuples[0], puweights, hrange)
    generate_files(arguments, nthreads)
    #for n in tqdm(ntuples):
    #    apply_corrections(n, x, mz_mc, mz_dt, pt_sf, mz_res_mc, mz_res_dt)
