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


def base_filename(path):
    return path.split("/")[-3]


def job_wrapper(args):
    # print(args)
    return friend_producer(*args)

def remove_zero_event_file(rfile):
    f = ROOT.TFile(rfile, "read")
    is_non_zero = f.GetListOfKeys().Contains("ntuple")

    if not is_non_zero:
        output_path_zero = rfile.replace("ntuples", "ntuples_zero_event")
        assert (not os.path.exists(output_path_zero))
        if not os.path.exists(os.path.dirname(output_path_zero)):
            os.makedirs(os.path.dirname(output_path_zero), exist_ok=False)
        shutil.move(rfile, output_path_zero)

    return is_non_zero

def friend_producer(rfile, dataset_proc):
    is_non_zero = remove_zero_event_file(rfile)
    if not is_non_zero:
        return

    output_path = rfile.replace("ntuples", "friends/xsec")

    if os.path.exists(output_path):
        print(f"friend_producer: {output_path} exists -> skip")
        return

    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=False)

    rdf = ROOT.RDataFrame("ntuple", rfile)

    if "Run20" in rfile:
        rdf = rdf.Define(
            "numberGeneratedEventsWeight",
            "(float)1.",
        )
        
        rdf = rdf.Define(
            "sumwWeight",
            "(float)1.",
        )

        rdf = rdf.Define(
            "sumwnormWeight",
            "(float)1.",
        )

        rdf = rdf.Define(
            "negFracWeight",
            "(float)0.",
        )

        rdf = rdf.Define(
            "crossSectionPerEventWeight",
            "(float)1.",
        )

        rdf = rdf.Define(
            "scale1fb_sumw",
            "(float)1.",
        )

        rdf = rdf.Define(
            "scale1fb_sumwnorm",
            "(float)1.",
        )

    else:
        rdf = rdf.Define("abs_genweight", "fabs(genweight)")
        weight = rdf.Mean("abs_genweight").GetValue()

        numberGeneratedEventsWeight = 1 / float(dataset_proc["nevents"])
        crossSectionPerEventWeight = float(dataset_proc["xsec"])
        sumwWeight = 1. / (float(dataset_proc["generator_weight"])*float(dataset_proc["nevents"])*weight)
        sumwnormWeight = weight * sumwWeight
        negFracWeight = 0.5 * (1-dataset_proc["generator_weight"])
        scale1fb_sumw = crossSectionPerEventWeight * sumwWeight * 1.e3
        scale1fb_sumwnorm = crossSectionPerEventWeight * sumwnormWeight * 1.e3

        rdf = rdf.Define(
            "numberGeneratedEventsWeight",
            "(float){ngw}".format(ngw=numberGeneratedEventsWeight),
        )
        
        rdf = rdf.Define(
            "sumwWeight",
            "(float){ngw}".format(ngw=sumwWeight),
        )

        rdf = rdf.Define(
            "sumwnormWeight",
            "(float){ngw}".format(ngw=sumwnormWeight),
        )

        rdf = rdf.Define(
            "negFracWeight",
            "(float){ngw}".format(ngw=negFracWeight),
        )

        rdf = rdf.Define(
            "crossSectionPerEventWeight",
            "(float){xsec}".format(xsec=crossSectionPerEventWeight),
        )

        rdf = rdf.Define(
            "scale1fb_sumw",
            "(float){ngw}".format(ngw=scale1fb_sumw),
        )

        rdf = rdf.Define(
            "scale1fb_sumwnorm",
            "(float){ngw}".format(ngw=scale1fb_sumwnorm),
        )

    rdf.Snapshot(
        "ntuple",
        output_path,
        [
            "numberGeneratedEventsWeight",
            "sumwWeight",
            "sumwnormWeight",
            "negFracWeight",
            "crossSectionPerEventWeight",
            "scale1fb_sumw",
            "scale1fb_sumwnorm",
        ],
    )


def generate_friend_trees(dataset, ntuples, nthreads):
    arguments = [(ntuple, dataset[base_filename(ntuple)]) for ntuple in ntuples]
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
    base_path = "/ceph/jdriesch/CROWN_samples/RerecoRun3V01/ntuples/*/*/*.root"
    dataset = yaml.load(open("datasets.yaml"), Loader=yaml.Loader)

    ntuples = glob.glob(base_path)
    ntuples_wo_data = ntuples.copy()
    nthreads = 64
    if nthreads > len(ntuples_wo_data):
        nthreads = len(ntuples_wo_data)
    print("# files: ", len(ntuples_wo_data))
    generate_friend_trees(dataset, ntuples_wo_data, nthreads)
