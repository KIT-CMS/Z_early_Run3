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
    
    rdf = rdf.Define('puweight', '0')
    quants = ["puweight", "puweightUp", "puweightDn"]
        
    for var in hrange.keys():
        bins = np.linspace(hrange[var][0], hrange[var][1], hrange[var][2]+1)

        rdf = rdf.Define(var+'weight','1')

        if not is_data:
            rdf = rdf.Redefine(var+'weight','1') # maybe better to initiate with 0
            for i in range(hrange[var][2]):
                rdf = rdf.Redefine(
                    var+'weight',
                    'double w;'\
                    f'if ({var} >= {bins[i]} && {var} < {bins[i+1]})'\
                    f'w = {wdict[var].GetBinContent(i)};'\
                    f'else w = {var}weight;'\
                    'return w;'
                )
            
        rdf = rdf.Redefine(
            'puweight', 
            f'puweight + 1./{len(hrange.keys())} * {var}weight'
            )
        quants += [var+'weight']

    # print(quants)
    mean = rdf.Mean("puweight").GetValue()
    rdf = rdf.Redefine("puweight", f"puweight/{mean}")
    rdf = rdf.Define("puweightUp", "puweight + abs(puweight-npvGoodweight)")
    mean_up = rdf.Mean("puweightUp").GetValue()
    rdf = rdf.Redefine("puweightUp", f"puweightUp/{mean_up}")
    rdf = rdf.Define("puweightDn", "puweight - abs(puweight-npvGoodweight)")
    mean_dn = rdf.Mean("puweightDn").GetValue()
    rdf = rdf.Redefine("puweightDn", f"puweightDn/{mean_dn}")
    
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
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    #ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    # ROOT.ROOT.EnableImplicitMT(16)
    ROOT.gROOT.SetBatch(True)


    # Load ntuples
    base_path = "/ceph/jdriesch/CROWN_samples/Run3V07/ntuples/2022/*/mm*/*.root"
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
    # print(arguments[0])
    # apply_corrections(ntuples[0], puweights, hrange)
    generate_files(arguments, nthreads)
    #for n in tqdm(ntuples):
    #    apply_corrections(n, x, mz_mc, mz_dt, pt_sf, mz_res_mc, mz_res_dt)
