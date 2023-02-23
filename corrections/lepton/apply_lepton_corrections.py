from ast import arg
import ROOT
import os
from tqdm import tqdm
import numpy as np
import array as a
import json
import glob
from multiprocessing import Pool, RLock

def calc_m(rdf, corr=""):
    rdf = rdf.Define(f"lv_1{corr}", f"ROOT::Math::PtEtaPhiMVector p(pt_1{corr}, eta_1, phi_1, mass_1); return p;")
    rdf = rdf.Define(f"lv_2{corr}", f"ROOT::Math::PtEtaPhiMVector p(pt_2{corr}, eta_2, phi_2, mass_2); return p;")
    
    rdf = rdf.Define(f"dilep{corr}", f"lv_1{corr} + lv_2{corr}")
    rdf = rdf.Define(f"pt_vis{corr}", f"dilep{corr}.Pt()")
    rdf = rdf.Define(f"m_vis{corr}", f"dilep{corr}.M()")
    rdf = rdf.Define(f"rap_vis{corr}", f"dilep{corr}.Rapidity()")
    return rdf


ROOT.gInterpreter.Declare("""
        float gaus(){
            return gRandom->Gaus(0,1);
        }
        """)

def job_wrapper(args):
    return apply_corrections(*args)

def apply_corrections(f, x, mz_mc, mz_dt, pt_sf, mz_res_mc, mz_res_dt):
    output_path = f.replace("ntuples_xsec_sf", "ntuples_xsec_sf_scaleres")
    
    rdf = ROOT.RDataFrame('ntuple', f)

    # extract original columns for later snapshot
    original_cols = [str(col) for col in rdf.GetColumnNames()]
    
    # check if data
    is_data = (rdf.Sum("is_data").GetValue()>0)
    
    # check if muon or electron event and set bins accordingly
    is_muonic = (('/mm/' in f) or ('/mmet/' in f))
    pt = [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]
    if is_muonic:
        eta = [-2.4, -1.2, 0., 1.2, 2.4]
    else:
        eta = [-2.5, -1.49, 0., 1.49, 2.5]
    
    # check if dimuon event and define corrected lepton quantities accordingly
    is_dilepton = ('/mm/' in f)
    rdf = rdf.Define("pt_1_corr", "pt_1")
    rdf = rdf.Define("pt_1_corr_up", "pt_1")
    rdf = rdf.Define("pt_1_corr_dn", "pt_1")
    if is_dilepton:
        rdf = rdf.Define("pt_2_corr", "pt_2")
        rdf = rdf.Define("pt_2_corr_up", "pt_2")
        rdf = rdf.Define("pt_2_corr_dn", "pt_2")
    
    for j in range(len(eta)-1): 
        for k in range(len(pt)-1): 
            eta_l, eta_r = eta[j], eta[j+1]
            pt_l, pt_r = pt[k], pt[k+1]

            # filter for a muon in event in corresponding bins
            filter1a = f"(eta_1 > {eta_l} && eta_1 < {eta_r})"
            filter1b = f"(pt_1 > {pt_l} && pt_1 < {pt_r})"
            # set resolution ratio to one if data resolution is larger
            if mz_res_mc[j][k] > mz_res_dt[j][k]:
                res_sf = 1
            else:
                res_sf = mz_res_dt[j][k] / mz_res_mc[j][k]
            # define momentum scale and resolution smearing
            mom_scale = f"(91.1876 + {mz_mc[j][k]})/(91.1876 + {mz_dt[j][k]})"
            res_smear = f"(1+{x[j][k]*mz_res_mc[j][k]}/91.1876*sqrt({res_sf**2}-1)*(float)(gaus()))"
            # in data events: scale momentum to match mass peak position in mc
            if is_data:
                rdf = rdf.Redefine(
                    "pt_1_corr", 
                    "double p;"\
                        f"if ({filter1a} && {filter1b}) p=pt_1 * {mom_scale};"\
                        "else p=pt_1_corr;"\
                        "return p;"
                    )
            # in mc events: smear pt to match mass resolution in data
            else:
                rdf = rdf.Redefine(
                    "pt_1_corr", 
                    "double p;"\
                        f"if ({filter1a} && {filter1b}) p=pt_1 * {res_smear};"\
                        "else p=pt_1_corr;"\
                        "return p;"
                    )
                rdf = rdf.Redefine(
                    "pt_1_corr_up",
                    "double p;"\
                        f"if ({filter1a} && {filter1b})"\
                        f"p=pt_1_corr * (1 + 0.5*abs({res_smear}-1) + 0.5*abs(1./({mom_scale})-1));"\
                        "else p=pt_1_corr_up;"\
                        "return p;"
                    )
                rdf = rdf.Redefine(
                    "pt_1_corr_dn",
                    "double p;"\
                        f"if ({filter1a} && {filter1b})"\
                        f"p=pt_1_corr * (1 - 0.5*abs({res_smear}-1) - 0.5*abs(1./({mom_scale})-1));"\
                        "else p=pt_1_corr_dn;"\
                        "return p;"
                    )


            # for z events: also correct second lepton
            if is_dilepton:
                filter2a = f"(eta_2 > {eta_l} && eta_2 < {eta_r})"
                filter2b = f"(pt_2 > {pt_l} && pt_2 < {pt_r})"
                
                if is_data:
                    rdf = rdf.Redefine(
                        "pt_2_corr", 
                        "double p;"\
                            f"if ({filter2a} && {filter2b}) p=pt_2 * {mom_scale};"\
                            "else p=pt_2_corr;"\
                            "return p;"
                        )
                else:
                    rdf = rdf.Redefine(
                        "pt_2_corr", 
                        "double p;"\
                            f"if ({filter2a} && {filter2b}) p=pt_2 * {res_smear};"\
                            "else p=pt_2_corr;"\
                            "return p;"
                        )
                    rdf = rdf.Redefine(
                        "pt_2_corr_up",
                        "double p;"\
                            f"if ({filter2a} && {filter2b})"\
                            f"p=pt_2_corr * (1 + 0.5*abs({res_smear}-1) + 0.5*abs(1./({mom_scale})-1));"\
                            "else p=pt_2_corr_up;"\
                            "return p;"
                        )
                    rdf = rdf.Redefine(
                        "pt_2_corr_dn",
                        "double p;"\
                            f"if ({filter2a} && {filter2b})"\
                            f"p=pt_2_corr * (1 - 0.5*abs({res_smear}-1) - 0.5*abs(1./({mom_scale})-1));"\
                            "else p=pt_2_corr_dn;"\
                            "return p;"
                        )
   # set variation to corrected values in data 
    if is_data:
        rdf = rdf.Redefine("pt_1_corr_up", "pt_1_corr")
        rdf = rdf.Redefine("pt_1_corr_dn", "pt_1_corr")
        if is_dilepton:
            rdf = rdf.Redefine("pt_2_corr_up", "pt_2_corr")
            rdf = rdf.Redefine("pt_2_corr_dn", "pt_2_corr")
                    
    # calculate corrected visible mass for dilepton events
    if is_dilepton:
        rdf = calc_m(rdf, '_corr')
        rdf = calc_m(rdf, '_corr_up')
        rdf = calc_m(rdf, '_corr_dn')
        quants =  [
                "pt_1_corr", "pt_1_corr_up", "pt_1_corr_dn",
                "pt_2_corr", "pt_2_corr_up", "pt_2_corr_dn",
                "m_vis_corr", "m_vis_corr_up", "m_vis_corr_dn",
                "pt_vis_corr", "pt_vis_corr_up", "pt_vis_corr_dn",
                "rap_vis_corr",
                ]
    else:
        rdf = rdf.Define("lv_1", "ROOT::Math::PtEtaPhiMVector p(pt_1_corr, eta_1, phi_1, mass_1); return p")
        rdf = rdf.Define("rap_vis_corr", "lv_1.Rapidity()")
        quants =  ["pt_1_corr", "pt_1_corr_up", "pt_1_corr_dn"]

    # correct and add recoil variables
    if not ("WtoLNu" in output_path):
        bosonphi = "phi_vis_c"
        bosonpt = "pt_vis_c"
        bosonrap = "rap_vis_corr"
    else:
        bosonphi = "genbosonphi"
        bosonpt = "genbosonpt"
        bosonrap = "genbosonrapidity"

    if is_dilepton:
        rdf = rdf.Define("pt_vis_c_x", "pt_1_corr*cos(phi_1) + pt_2_corr*cos(phi_2)")
        rdf = rdf.Define("pt_vis_c_y", "pt_1_corr*sin(phi_1) + pt_2_corr*sin(phi_2)")
        rdf = rdf.Define("pt_vis_c", "sqrt(pt_vis_c_x*pt_vis_c_x + pt_vis_c_y*pt_vis_c_y)")
        rdf = rdf.Define("phi_vis_c", "atan2(pt_vis_c_y, pt_vis_c_x)")
    else:
        rdf = rdf.Define("pt_vis_c", "pt_1_corr")
        rdf = rdf.Define("phi_vis_c", "phi_1")

    rdf = rdf.Define("bosonpt", f"{bosonpt}")
    rdf = rdf.Define("bosonphi", f"{bosonphi}")
    rdf = rdf.Define("bosonrap", f"{bosonrap}")

    rdf = rdf.Define("uPx", "met_uncorrected*cos(metphi_uncorrected) + pt_vis_c*cos(phi_vis_c)")
    rdf = rdf.Define("uPy", "met_uncorrected*sin(metphi_uncorrected) + pt_vis_c*sin(phi_vis_c)")

    rdf = rdf.Define("uP1_uncorrected", f"- (uPx*cos({bosonphi}) + uPy*sin({bosonphi}))")
    rdf = rdf.Define("uP2_uncorrected", f"uPx*sin({bosonphi}) - uPy*cos({bosonphi})")

    rdf = rdf.Define("pfuPx", "pfmet_uncorrected*cos(pfmetphi_uncorrected) + pt_vis_c*cos(phi_vis_c)")
    rdf = rdf.Define("pfuPy", "pfmet_uncorrected*sin(pfmetphi_uncorrected) + pt_vis_c*sin(phi_vis_c)")

    rdf = rdf.Define("pfuP1_uncorrected", f"- (pfuPx*cos({bosonphi}) + pfuPy*sin({bosonphi}))")
    rdf = rdf.Define("pfuP2_uncorrected", f"pfuPx*sin({bosonphi}) - pfuPy*cos({bosonphi})")

    met_cols = [
        "uP1_uncorrected", "uP2_uncorrected", "pfuP1_uncorrected", "pfuP2_uncorrected", 
        "pt_vis_c", "phi_vis_c", "bosonpt", "bosonphi", "bosonrap"]

    outdir = output_path.replace(output_path.split('/')[-1],"")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    #print(original_cols + quants + met_cols)
    rdf.Snapshot("ntuple", output_path, original_cols + quants + met_cols)

    #print("Great success!")
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


    # base_path = "/work/jdriesch/earlyrun3/samples/Run3V04/ntuples_xsec_sf_EraC/20*/*/*/*.root"
    base_path = "/storage/9/jdriesch/earlyrun3/samples/Run3V06/ntuples_xsec_sf_EraC/20*/*/*/*.root"
    ntuples = glob.glob(base_path)
    # Load correction files
    x = np.loadtxt('correction_files/Run3/mm/res_sf_extra.txt')
    mz_mc = np.loadtxt('correction_files/Run3/mm/mz_mc.txt')
    mz_dt = np.loadtxt('correction_files/Run3/mm/mz_dt.txt')
    pt_sf = (91.1876+mz_mc) / (91.1876+mz_dt)
    mz_res_mc = np.loadtxt('correction_files/Run3/mm/resmz_mc.txt')
    mz_res_dt = np.loadtxt('correction_files/Run3/mm/resmz_dt.txt') * pt_sf
    
    nthreads = 16
    arguments = [(ntuple, x, mz_mc, mz_dt, pt_sf, mz_res_mc, mz_res_dt) for ntuple in ntuples]
    generate_files(arguments, nthreads)
    #for n in tqdm(ntuples):
    #    apply_corrections(n, x, mz_mc, mz_dt, pt_sf, mz_res_mc, mz_res_dt)
