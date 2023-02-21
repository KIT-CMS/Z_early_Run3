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

def apply_corrections(f, puweights, hrange):
    output_path = f.replace("ntuples_xsec_sf_scaleres", "ntuples_xsec_sf_scaleres_pu")
    
    rdf = ROOT.RDataFrame('ntuple', f)

    # extract original columns for later snapshot
    original_cols = [str(col) for col in rdf.GetColumnNames()]
    
    # check if data
    is_data = (rdf.Sum("is_data").GetValue()>0) 
    
    rdf = rdf.Define('puweight', '0')
    quants = original_cols + ["puweight"]
        
    for var in hrange.keys():
        bins = np.linspace(hrange[var][0], hrange[var][1], hrange[var][2]+1)

        rdf = rdf.Define(var+'weight','1')
        
        if not is_data:
            for i in range(hrange[var][2]):
                rdf = rdf.Redefine(
                    var+'weight',
                    'double w;'\
                    f'if ({var} > {bins[i]} && {var} < {bins[i+1]})'\
                    f'w = {puweights[var][i]};'\
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
    ROOT.ROOT.EnableImplicitMT(16)
    #ROOT.gROOT.SetBatch(True)


    # Load ntuples
    base_path = "/storage/9/jdriesch/earlyrun3/samples/Run3V06/ntuples_xsec_sf_scaleres_EraC/20*/*/*/*.root"
    ntuples = glob.glob(base_path)

    # Load correction files
    jsonfile = 'puweights.json'
    with open(jsonfile) as f:
        puweights = json.load(f)

    hrange = {
        'npvGood': [0, 60, 60],
        'rhoFastjetCentralChargedPileUp': [0, 50, 50],
        'rhoFastjetCentralCalo': [0, 35, 35]
    }    
    nthreads = 16
    arguments = [(ntuple, puweights, hrange) for ntuple in ntuples]
    # print(arguments[0])
    # apply_corrections(ntuples[0], puweights, hrange)
    generate_files(arguments, nthreads)
    #for n in tqdm(ntuples):
    #    apply_corrections(n, x, mz_mc, mz_dt, pt_sf, mz_res_mc, mz_res_dt)
