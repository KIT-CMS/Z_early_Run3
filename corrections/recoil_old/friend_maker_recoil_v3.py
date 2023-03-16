#from sys import base_prefix
import ROOT
import argparse
import os
from tqdm import tqdm
import numpy as np
import glob
from multiprocessing import Pool, current_process, RLock


def parse_args():
    parser = argparse.ArgumentParser(description="produce friend of input ntuple with corrected met")   
    parser.add_argument('-I', '--inpath', default='/ceph/moh/CROWN_samples/Run3V01/ntuples/2022/*/mm*/*.root', help='path to input samples')
    parser.add_argument('--overwrite', action='store_true', default=True)
    parser.add_argument('--test', action='store_true', default=False)   
    args = parser.parse_args()
    return args


def load_fit_results(pathName, nBins):
    # loop over parallel and perpendicular component of hadronic recoil and load fit results

    workspace = []

    for i in ["1", "2"]:
        file = ROOT.TFile(pathName+"pdfsU" + i + ".root")
        workspace.append(file.Get("pdfsU"+i))

        for j in range(nBins):
            pdf = workspace[-1].pdf("sig_" + str(j))
            myX = workspace[-1].var("u_" + str(j))
            cdf = pdf.createCdf(myX)
            workspace[-1].Import(cdf, ROOT.RooFit.RecycleConflictNodes() ,ROOT.RooFit.Silence())    

    # TODO: implement statistical uncertainty calculation?
    return workspace


def invertCdf(uP, mc_cdf, mc_x, data_cdf, data_x):
    # calculate p-value of mc distribution
    if (uP < mc_x.getMin()) or (uP > mc_x.getMax()): 
        return uP
    mc_x.setVal(uP)
    data_x.setVal(uP)
    data_cdf.getVal()
    # calculate inverse value of mc p-value
    pVal = data_cdf.findRoot(data_x, data_x.getMin(), data_x.getMax(), mc_cdf.getVal())
    data_x.setVal(pVal)
    return pVal


def prep(inpath, met="met"):    
    # load mc to correct
    chain = ROOT.TChain("ntuple")
    chain.Add(inpath)
    friend = ROOT.TChain("ntuple")
    frpath = inpath.replace("ntuples", "friends/lep_corr_v3")
    # print("prep inpath:", inpath)
    # print("prep frpath:", frpath)
    friend.Add(frpath)
    chain.AddFriend(friend)
    rdf = ROOT.RDataFrame(chain)
    # rdf = ROOT.RDataFrame("ntuple", inpath)

    # check how many leptons in final state and adjust total leptonic momentum accordingly
    channel = inpath.split("/")[-2]

    # apply rochester corrections if muon involved (mm or mmet)
    # if "mm" in channel:
    #     cor = "_corr"
    # else:
    #     cor = ""
    cor = "_corr"

    # if not ("DYtoLL" in inpath or "WtoLNu" in inpath):
    if not ("WtoLNu" in inpath):
        bosonphi = "phi_vis_c"
        bosonpt = "pt_vis_c"
    else:
        bosonphi = "genbosonphi"
        bosonpt = "genbosonpt"

    if channel == "ee" or channel == "mm":
        rdf = rdf.Define("pt_vis_c_x", "pt_1"+cor+"*cos(phi_1) + pt_2"+cor+"*cos(phi_2)")
        rdf = rdf.Define("pt_vis_c_y", "pt_1"+cor+"*sin(phi_1) + pt_2"+cor+"*sin(phi_2)")
        rdf = rdf.Define("pt_vis_c", "sqrt(pt_vis_c_x*pt_vis_c_x + pt_vis_c_y*pt_vis_c_y)")
        rdf = rdf.Define("phi_vis_c", "atan2(pt_vis_c_y, pt_vis_c_x)")
    else:
        rdf = rdf.Define("pt_vis_c", "pt_1"+cor+"")
        rdf = rdf.Define("phi_vis_c", "phi_1")

    rdf = rdf.Define("uPx", met+"_uncorrected*cos("+met+"phi_uncorrected) + pt_vis_c*cos(phi_vis_c)") # sign is wrong but consistently -> met is correct
    rdf = rdf.Define("uPy", met+"_uncorrected*sin("+met+"phi_uncorrected) + pt_vis_c*sin(phi_vis_c)")

    rdf = rdf.Define("uP1_uncorrected", "- (uPx*cos("+bosonphi+") + uPy*sin("+bosonphi+"))") # maybe use reco-phi instead
    rdf = rdf.Define("uP2_uncorrected", "uPx*sin("+bosonphi+") - uPy*cos("+bosonphi+")")

    return rdf, chain, friend


def calculateMET(rdf, corrected, data=False, met="met"):
    rdf = rdf.Define(met+"Px"+corrected, "-pt_vis_c*cos(phi_vis_c) - uP1"+corrected+"*cos(genbosonphi) + uP2"+corrected+"*sin(genbosonphi)") # TODO check if name already used
    rdf = rdf.Define(met+"Py"+corrected, "-pt_vis_c*sin(phi_vis_c) - uP1"+corrected+"*sin(genbosonphi) - uP2"+corrected+"*cos(genbosonphi)")

    rdf = rdf.Define(met+corrected, "sqrt("+met+"Px"+corrected+"*"+met+"Px"+corrected+" + "+met+"Py"+corrected+"*"+met+"Py"+corrected+")")
    rdf = rdf.Define(met+"phi"+corrected, "atan2("+met+"Py"+corrected+", "+met+"Px"+corrected+")")

    return rdf


def run(input_dict):
    infile = input_dict["ntuple"]
    args = input_dict["args"]
    ws_dict = input_dict["ws_dict"]

    version = infile.split('/')[-6]
    process = infile.split('/')[-3]
    year = infile.split('/')[-4]
    channel = infile.split('/')[-2]

    zPtBinEdges = [0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000]
    nBins = len(zPtBinEdges)-1

    if channel == 'mmet':
        dataWorkspace = load_fit_results(ws_dict["data_mm"], nBins)
        mcWorkspace_pos = load_fit_results(ws_dict["Zmm"], nBins)
        mcWorkspace_neg = load_fit_results(ws_dict["Zmm"], nBins)           
        mcWorkspace_pos_tocorr = load_fit_results(ws_dict["Wmpos"], nBins)
        mcWorkspace_neg_tocorr = load_fit_results(ws_dict["Wmneg"], nBins)   

    elif channel == 'emet':
        dataWorkspace = load_fit_results(ws_dict["data_ee"], nBins)
        mcWorkspace_pos = load_fit_results(ws_dict["Zee"], nBins)
        mcWorkspace_neg = load_fit_results(ws_dict["Zee"], nBins)
        mcWorkspace_pos_tocorr = load_fit_results(ws_dict["Wepos"], nBins)
        mcWorkspace_neg_tocorr = load_fit_results(ws_dict["Weneg"], nBins)

    elif channel == 'mm':
        dataWorkspace = load_fit_results(ws_dict["data_mm"], nBins)
        mcWorkspace = load_fit_results(ws_dict["Zmm"], nBins)
        mcWorkspace_tocorr = load_fit_results(ws_dict["Zmm"], nBins)

    elif channel == 'ee':
        dataWorkspace = load_fit_results(ws_dict["data_ee"], nBins)
        mcWorkspace = load_fit_results(ws_dict["Zee"], nBins)
        mcWorkspace_tocorr = load_fit_results(ws_dict["Zee"], nBins)

    else:
        print('channel cannot be matched to mm, ee, mmet, emet')
        return

    # make sure that variables are named consistently (for Run2 met is already corrected)
    if "18" in year:
        corrected = "_corrected"
    else:
        corrected = "_corr"
    met = 'met'
    if 'pfmet_' in ws_dict["data_mm"]:
        met = 'pfmet'
    
    if args.test:
        outdir = '/ceph/moh/CROWN_samples/test/'
    else:
        outdir = infile.replace("/ntuples/", "/friends/"+met+"_corr_with_lep_corr_v3/")

    outfile = outdir
    if args.test:
        outfile = outdir + infile.split("/")[-1]

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        # print("Created new directory: ", outdir)
        
    rdf, _chain, _friend = prep(infile, met)

    # for other processes than DY: just use uncorrected values
    if not ("DYtoLL" in process or "WtoLNu" in process):
        rdf = rdf.Define("uP1"+corrected, "uP1_uncorrected")
        rdf = rdf.Define("uP2"+corrected, "uP2_uncorrected")
        rdf = rdf.Define(met+corrected, 'double('+met+'_uncorrected)')
        rdf = rdf.Define(met+"phi"+corrected, 'double('+met+'phi_uncorrected)')

        if not os.path.exists(outfile) or args.overwrite:
            rdf.Snapshot("ntuple", outfile, ["uP1_uncorrected", "uP2_uncorrected", "uP1"+corrected, "uP2"+corrected, met+corrected, met+"phi"+corrected, "pt_vis_c", "phi_vis_c"])

        return

    rdf = rdf.Define("uP1"+corrected, "uP1_uncorrected")
    rdf = rdf.Define("uP2"+corrected, "uP2_uncorrected")

    columns = [
        "uP1"+corrected,
        "uP2"+corrected,
        "uP1_uncorrected",
        "uP2_uncorrected",
        "pt_vis_c",
        "phi_vis_c",
        "genbosonpt",
        "genbosonphi",
        'q_1'
    ]
    fromrdf = rdf.AsNumpy(columns = columns)
    tordf = {
        "uP1"+corrected: fromrdf["uP1"+corrected], 
        "uP2"+corrected: fromrdf["uP2"+corrected],             
        "uP1_uncorrected": fromrdf["uP1_uncorrected"], 
        "uP2_uncorrected": fromrdf["uP2_uncorrected"], 
        "pt_vis_c": fromrdf["pt_vis_c"], 
        "phi_vis_c": fromrdf["phi_vis_c"], 
        "genbosonpt": fromrdf["genbosonpt"], 
        "genbosonphi": fromrdf["genbosonphi"]
    }

    # for j in tqdm(range(len(fromrdf["genbosonpt"]))):
    for j in range(len(fromrdf["genbosonpt"])):
        # Get pT bin of genboson
        kBin = nBins-1
        for k in range(nBins):
            if ((channel == 'mm' or channel == 'ee') and fromrdf["pt_vis_c"][j] < zPtBinEdges[k]) or ((channel == 'mmet' or channel == 'emet') and fromrdf["genbosonpt"][j] < zPtBinEdges[k]): # TODO check if right order
                kBin = k-1
                break
        sBin = str(kBin)

        # MET corrections for parallel and perpendicular recoil component
        for k in [0,1]:
            """
            data_pdf   = dataWorkspace[k].pdf('sig_'+sBin)
            mc_pdf     = mcWorkspace[k].pdf('sig_'+sBin)
            toCorr_pdf = toCorrWorkspace[k].pdf('sig_'+sBin)
            """
            
            data_cdf   = dataWorkspace[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
            data_xuP   = dataWorkspace[k].var('u_'+sBin)

            if ('met' in channel) and (fromrdf["q_1"][j] == 1):
                mc_cdf     = mcWorkspace_pos[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf = mcWorkspace_pos_tocorr[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')    
                mc_xuP     = mcWorkspace_pos[k].var('u_'+sBin)
                toCorr_xuP = mcWorkspace_pos_tocorr[k].var('u_'+sBin)

            elif ('met' in channel) and (fromrdf["q_1"][j] == -1):                   
                mc_cdf     = mcWorkspace_neg[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf = mcWorkspace_neg_tocorr[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                mc_xuP     = mcWorkspace_neg[k].var('u_'+sBin)
                toCorr_xuP = mcWorkspace_neg_tocorr[k].var('u_'+sBin)

            else:
                mc_cdf     = mcWorkspace[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf = mcWorkspace_tocorr[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                mc_xuP     = mcWorkspace[k].var('u_'+sBin)
                toCorr_xuP = mcWorkspace_tocorr[k].var('u_'+sBin)


            mc_uPValZlike = invertCdf(fromrdf["uP"+str(k+1)+"_uncorrected"][j], toCorr_cdf, toCorr_xuP, mc_cdf, mc_xuP) # TODO check order of matrix
            dt_uPValZlike = invertCdf(fromrdf["uP"+str(k+1)+"_uncorrected"][j], toCorr_cdf, toCorr_xuP, data_cdf, data_xuP)
            # data_uPValZlike = invertCdf(mc_uPValZlike, mc_cdf, mc_xuP, data_cdf, data_xuP)

            tordf["uP"+str(k+1)+corrected][j] += dt_uPValZlike - mc_uPValZlike

    rdf_tosave = ROOT.RDF.MakeNumpyDataFrame(tordf)
    rdf_tosave = calculateMET(rdf_tosave, corrected, met=met)     

    if not os.path.exists(outfile) or args.overwrite:
        rdf_tosave.Snapshot("ntuple", outfile, ["uP1_uncorrected", "uP2_uncorrected", "uP1"+corrected, "uP2"+corrected, met+corrected, met+"phi"+corrected, "pt_vis_c", "phi_vis_c"])

    print(">>> "+infile+" success")


if __name__=='__main__':
    args = parse_args()

    ws_dicts = [
        {
            "data_mm":  "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V01_outputs/met_uncorrected_data_triple_muon_sigAndBck/", 
            "Zmm":      "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V01_outputs/met_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V01_outputs/met_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V01_outputs/met_uncorrected_WposMC_triple_muon_sigOnly/",
        },
        {
            "data_mm":  "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V01_outputs/pfmet_uncorrected_data_triple_muon_sigAndBck/", 
            "Zmm":      "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V01_outputs/pfmet_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V01_outputs/pfmet_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V01_outputs/pfmet_uncorrected_WposMC_triple_muon_sigOnly/",
        }
    ]

    ntuples = glob.glob(args.inpath)
    nthreads = 64

    arguments = [{"ntuple": ntuple, "args": args, "ws_dict": ws_dict} for ntuple in ntuples for ws_dict in ws_dicts]
    pool = Pool(nthreads, initargs=(RLock(),), initializer=tqdm.set_lock)
    for _ in tqdm(
        pool.imap_unordered(run, arguments),  # imap_unordered
        total=len(arguments),
        desc="Total progess",
        # position=nthreads + 1,
        dynamic_ncols=True,
        leave=True,
    ):
        pass
