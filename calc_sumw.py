from operator import ne
import ROOT
import argparse
import yaml
import os, sys
import time
import glob
from tqdm import tqdm
from multiprocessing import Pool, current_process, RLock

import uproot
import pandas as pd

def base_filename(path):
    return path.split("/")[-1]

def job_wrapper(args):
    # print(args)
    return calc(*args)

def calc(path, dataset_proc):
    df = uproot.concatenate(f'{path}/mm/*.root:ntuple', num_workers = 8, library='pd')
    nevents = len(df.index)

    # assert nevents == dataset_proc["nevents"]

    sumw = df['genweight'].astype('float64').sum()

    df.loc[:,'genweight_sign'] = df['genweight'].astype('float64') / abs(df['genweight'].astype('float64'))
    sumwnorm = df['genweight_sign'].astype('float64').sum()
    negfrac = (0.5*(1. - (sumwnorm/nevents)))
    generator_weight = (1. - 2.*negfrac)

    outstr = f"""
{base_filename(path)}
  nevents: {nevents:.10f}
  sumw: {sumw:.10f}
  sumwnorm: {sumwnorm:.10f}
  negfrac: {negfrac:.10f}
  generator_weight: {generator_weight:.10f}
"""

    print(outstr)
    sys.stdout.flush()

def calc_sumw(dataset, ntuples, nthreads):
    arguments = [(ntuple, dataset[base_filename(ntuple)]) for ntuple in ntuples]
    pool = Pool(nthreads, initargs=(RLock(),), initializer=tqdm.set_lock)
    for _ in pool.imap_unordered(job_wrapper, arguments):
        pass

if __name__ == "__main__":
    base_path = "/ceph/jdriesch/CROWN_samples/Run3V07_sumw/ntuples/2022/*"
    dataset = yaml.load(open("datasets.yaml"), Loader=yaml.Loader)

    ntuples = glob.glob(base_path)
    ntuples_wo_data = ntuples.copy()
    for ntuple in ntuples:
        if "Run20" in ntuple:
            ntuples_wo_data.remove(str(ntuple))
    nthreads = 16
    if nthreads > len(ntuples_wo_data):
        nthreads = len(ntuples_wo_data)
    calc_sumw(dataset, ntuples_wo_data, nthreads)
