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
    parser.add_argument('-I', '--inpath', default='/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC_lep_corr_01_x0p60/2022/*/mm*/*.root', help='path to input samples')
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
            workspace[-1].Import(cdf, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())    

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


def prep(inpath):    
    # load mc to correct
    chain = ROOT.TChain("ntuple")
    chain.Add(inpath)
    # friend = ROOT.TChain("ntuple")
    # frpath = inpath.replace("ntuples", "friends/lep_corr_v01_x0p7")
    # friend.Add(frpath)
    # chain.AddFriend(friend)
    rdf = ROOT.RDataFrame(chain)
    # rdf = ROOT.RDataFrame("ntuple", inpath)

    # # check how many leptons in final state and adjust total leptonic momentum accordingly
    # channel = inpath.split("/")[-2]

    # # apply rochester corrections if muon involved (mm or mmet)
    # # if "mm" in channel:
    # #     cor = "_corr"
    # # else:
    # #     cor = ""
    # cor = "_corr"

    # # if not ("DYtoLL" in inpath or "WtoLNu" in inpath):
    # if not ("WtoLNu" in inpath):
    #     bosonphi = "phi_vis_c"
    #     bosonpt = "pt_vis_c"
    # else:
    #     bosonphi = "genbosonphi"
    #     bosonpt = "genbosonpt"

    # if channel == "ee" or channel == "mm":
    #     rdf = rdf.Define("pt_vis_c_x", "pt_1"+cor+"*cos(phi_1) + pt_2"+cor+"*cos(phi_2)")
    #     rdf = rdf.Define("pt_vis_c_y", "pt_1"+cor+"*sin(phi_1) + pt_2"+cor+"*sin(phi_2)")
    #     rdf = rdf.Define("pt_vis_c", "sqrt(pt_vis_c_x*pt_vis_c_x + pt_vis_c_y*pt_vis_c_y)")
    #     rdf = rdf.Define("phi_vis_c", "atan2(pt_vis_c_y, pt_vis_c_x)")
    # else:
    #     rdf = rdf.Define("pt_vis_c", "pt_1"+cor+"")
    #     rdf = rdf.Define("phi_vis_c", "phi_1")

    # rdf = rdf.Define("uPx", met+"_uncorrected*cos("+met+"phi_uncorrected) + pt_vis_c*cos(phi_vis_c)") # sign is wrong but consistently -> met is correct
    # rdf = rdf.Define("uPy", met+"_uncorrected*sin("+met+"phi_uncorrected) + pt_vis_c*sin(phi_vis_c)")

    # rdf = rdf.Define("uP1_uncorrected", "- (uPx*cos("+bosonphi+") + uPy*sin("+bosonphi+"))") # maybe use reco-phi instead
    # rdf = rdf.Define("uP2_uncorrected", "uPx*sin("+bosonphi+") - uPy*cos("+bosonphi+")")

    return rdf, chain  # , friend


def calculateMET(rdf, corrected, data=False, mettype="", syst_postfix = ""):
    rdf = rdf.Define(mettype+"metPx"+corrected+syst_postfix, "-pt_vis_c*cos(phi_vis_c) - "+mettype+"uP1"+corrected+syst_postfix+"*cos(genbosonphi) + "+mettype+"uP2"+corrected+syst_postfix+"*sin(genbosonphi)")
    rdf = rdf.Define(mettype+"metPy"+corrected+syst_postfix, "-pt_vis_c*sin(phi_vis_c) - "+mettype+"uP1"+corrected+syst_postfix+"*sin(genbosonphi) - "+mettype+"uP2"+corrected+syst_postfix+"*cos(genbosonphi)")

    rdf = rdf.Define(mettype+"met"+corrected+syst_postfix, "sqrt("+mettype+"metPx"+corrected+syst_postfix+"*"+mettype+"metPx"+corrected+syst_postfix+" + "+mettype+"metPy"+corrected+syst_postfix+"*"+mettype+"metPy"+corrected+syst_postfix+")")
    rdf = rdf.Define(mettype+"metphi"+corrected+syst_postfix, "atan2("+mettype+"metPy"+corrected+syst_postfix+", "+mettype+"metPx"+corrected+syst_postfix+")")

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
        mcWorkspace = load_fit_results(ws_dict["Zmm"], nBins)
        mcWorkspace_pos_tocorr = load_fit_results(ws_dict["Wmpos"], nBins)
        mcWorkspace_neg_tocorr = load_fit_results(ws_dict["Wmneg"], nBins)
        dataWorkspace_pf = load_fit_results(ws_dict["data_mm_pf"], nBins)
        mcWorkspace_pf = load_fit_results(ws_dict["Zmm_pf"], nBins)
        mcWorkspace_pos_tocorr_pf = load_fit_results(ws_dict["Wmpos_pf"], nBins)
        mcWorkspace_neg_tocorr_pf = load_fit_results(ws_dict["Wmneg_pf"], nBins)   

    elif channel == 'mm':
        dataWorkspace = load_fit_results(ws_dict["data_mm"], nBins)
        mcWorkspace = load_fit_results(ws_dict["Zmm"], nBins)
        mcWorkspace_tocorr = load_fit_results(ws_dict["Zmm"], nBins)
        dataWorkspace_pf = load_fit_results(ws_dict["data_mm_pf"], nBins)
        mcWorkspace_pf = load_fit_results(ws_dict["Zmm_pf"], nBins)
        mcWorkspace_tocorr_pf = load_fit_results(ws_dict["Zmm_pf"], nBins)

    # elif channel == 'emet':
    #     dataWorkspace = load_fit_results(ws_dict["data_ee"], nBins)
    #     mcWorkspace_pos = load_fit_results(ws_dict["Zee"], nBins)
    #     mcWorkspace_neg = load_fit_results(ws_dict["Zee"], nBins)
    #     mcWorkspace_pos_tocorr = load_fit_results(ws_dict["Wepos"], nBins)
    #     mcWorkspace_neg_tocorr = load_fit_results(ws_dict["Weneg"], nBins)

    # elif channel == 'ee':
    #     dataWorkspace = load_fit_results(ws_dict["data_ee"], nBins)
    #     mcWorkspace = load_fit_results(ws_dict["Zee"], nBins)
    #     mcWorkspace_tocorr = load_fit_results(ws_dict["Zee"], nBins)

    else:
        print('channel cannot be matched to mm, ee, mmet, emet')
        return

    # make sure that variables are named consistently (for Run2 met is already corrected)
    if "18" in year:
        corrected = "_corrected"
    else:
        corrected = "_corr"

    syst_postfix = ""
    if "_sigOnly" in ws_dict["data_mm"]:
        syst_postfix = "_sigOnly"
    elif "_double" in ws_dict["data_mm"]:
        syst_postfix = "_double"
    elif "_CB" in ws_dict["data_mm"]:
        syst_postfix = "_CB"

    if args.test:
        outdir = '/ceph/moh/CROWN_samples/test/'
    else:
        outdir = infile.replace("/ntuples_xsec_sf_EraC_lep_corr_01_x0p60/", "/friend_xsec_sf_EraC_lep_corr_01_x0p60_met_corr"+syst_postfix+"/")

    outfile = outdir
    if args.test:
        outfile = outdir + infile.split("/")[-1]

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        # print("Created new directory: ", outdir)
        
    rdf, _chain = prep(infile)
    original_cols = [str(col) for col in rdf.GetColumnNames()]

    # for other processes than DY: just use uncorrected values
    if not ("DYtoLL" in process or "WtoLNu" in process):
        rdf = rdf.Define("uP1"+corrected+syst_postfix, "uP1_uncorrected")
        rdf = rdf.Define("uP2"+corrected+syst_postfix, "uP2_uncorrected")
        rdf = rdf.Define("met"+corrected+syst_postfix, 'double(met_uncorrected)')
        rdf = rdf.Define("metphi"+corrected+syst_postfix, 'double(metphi_uncorrected)')

        rdf = rdf.Define("pfuP1"+corrected+syst_postfix, "pfuP1_uncorrected")
        rdf = rdf.Define("pfuP2"+corrected+syst_postfix, "pfuP2_uncorrected")
        rdf = rdf.Define("pfmet"+corrected+syst_postfix, 'double(pfmet_uncorrected)')
        rdf = rdf.Define("pfmetphi"+corrected+syst_postfix, 'double(pfmetphi_uncorrected)')

        if not os.path.exists(outfile) or args.overwrite:
            rdf.Snapshot("ntuple", outfile, [
                "uP1"+corrected+syst_postfix, "uP2"+corrected+syst_postfix, "met"+corrected+syst_postfix, "metphi"+corrected+syst_postfix,
                "pfuP1"+corrected+syst_postfix, "pfuP2"+corrected+syst_postfix, "pfmet"+corrected+syst_postfix, "pfmetphi"+corrected+syst_postfix
            ])
            # rdf.Snapshot("ntuple", outfile, list(set(original_cols + [
            #     "uP1"+corrected+syst_postfix, "uP2"+corrected+syst_postfix, "met"+corrected+syst_postfix, "metphi"+corrected+syst_postfix,
            #     "pfuP1"+corrected+syst_postfix, "pfuP2"+corrected+syst_postfix, "pfmet"+corrected+syst_postfix, "pfmetphi"+corrected+syst_postfix
            # ])))

        return

    rdf = rdf.Define("uP1"+corrected+syst_postfix, "uP1_uncorrected")
    rdf = rdf.Define("uP2"+corrected+syst_postfix, "uP2_uncorrected")
    rdf = rdf.Define("pfuP1"+corrected+syst_postfix, "pfuP1_uncorrected")
    rdf = rdf.Define("pfuP2"+corrected+syst_postfix, "pfuP2_uncorrected")

    columns = [
        "uP1"+corrected+syst_postfix,
        "uP2"+corrected+syst_postfix,
        "uP1_uncorrected",
        "uP2_uncorrected",
        "pfuP1"+corrected+syst_postfix,
        "pfuP2"+corrected+syst_postfix,
        "pfuP1_uncorrected",
        "pfuP2_uncorrected",
        "pt_vis_c",
        "phi_vis_c",
        "genbosonpt",
        "genbosonphi",
        'q_1'
    ]
    fromrdf = rdf.AsNumpy(columns = columns)
    tordf = {
        "uP1"+corrected+syst_postfix: fromrdf["uP1"+corrected+syst_postfix],
        "uP2"+corrected+syst_postfix: fromrdf["uP2"+corrected+syst_postfix],
        "uP1_uncorrected": fromrdf["uP1_uncorrected"],
        "uP2_uncorrected": fromrdf["uP2_uncorrected"],
        "pfuP1"+corrected+syst_postfix: fromrdf["pfuP1"+corrected+syst_postfix],
        "pfuP2"+corrected+syst_postfix: fromrdf["pfuP2"+corrected+syst_postfix],
        "pfuP1_uncorrected": fromrdf["pfuP1_uncorrected"],
        "pfuP2_uncorrected": fromrdf["pfuP2_uncorrected"],
        "pt_vis_c": fromrdf["pt_vis_c"],
        "phi_vis_c": fromrdf["phi_vis_c"],
        "genbosonpt": fromrdf["genbosonpt"],
        "genbosonphi": fromrdf["genbosonphi"]
    }

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

            # PUPPI
            data_cdf   = dataWorkspace[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
            data_xuP   = dataWorkspace[k].var('u_'+sBin)

            if ('met' in channel) and (fromrdf["q_1"][j] == 1):
                mc_cdf     = mcWorkspace[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf = mcWorkspace_pos_tocorr[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')    
                mc_xuP     = mcWorkspace[k].var('u_'+sBin)
                toCorr_xuP = mcWorkspace_pos_tocorr[k].var('u_'+sBin)

            elif ('met' in channel) and (fromrdf["q_1"][j] == -1):                   
                mc_cdf     = mcWorkspace[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf = mcWorkspace_neg_tocorr[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                mc_xuP     = mcWorkspace[k].var('u_'+sBin)
                toCorr_xuP = mcWorkspace_neg_tocorr[k].var('u_'+sBin)

            else:
                mc_cdf     = mcWorkspace[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf = mcWorkspace_tocorr[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                mc_xuP     = mcWorkspace[k].var('u_'+sBin)
                toCorr_xuP = mcWorkspace_tocorr[k].var('u_'+sBin)


            mc_uPValZlike = invertCdf(fromrdf["uP"+str(k+1)+"_uncorrected"][j], toCorr_cdf, toCorr_xuP, mc_cdf, mc_xuP) # TODO check order of matrix
            dt_uPValZlike = invertCdf(fromrdf["uP"+str(k+1)+"_uncorrected"][j], toCorr_cdf, toCorr_xuP, data_cdf, data_xuP)

            tordf["uP"+str(k+1)+corrected+syst_postfix][j] += dt_uPValZlike - mc_uPValZlike

            # PF
            data_cdf_pf   = dataWorkspace_pf[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
            data_xuP_pf   = dataWorkspace_pf[k].var('u_'+sBin)

            if ('met' in channel) and (fromrdf["q_1"][j] == 1):
                mc_cdf_pf     = mcWorkspace_pf[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf_pf = mcWorkspace_pos_tocorr_pf[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')    
                mc_xuP_pf     = mcWorkspace_pf[k].var('u_'+sBin)
                toCorr_xuP_pf = mcWorkspace_pos_tocorr_pf[k].var('u_'+sBin)

            elif ('met' in channel) and (fromrdf["q_1"][j] == -1):                   
                mc_cdf_pf     = mcWorkspace_pf[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf_pf = mcWorkspace_neg_tocorr_pf[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                mc_xuP_pf     = mcWorkspace_pf[k].var('u_'+sBin)
                toCorr_xuP_pf = mcWorkspace_neg_tocorr_pf[k].var('u_'+sBin)

            else:
                mc_cdf_pf     = mcWorkspace_pf[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                toCorr_cdf_pf = mcWorkspace_tocorr_pf[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                mc_xuP_pf     = mcWorkspace_pf[k].var('u_'+sBin)
                toCorr_xuP_pf = mcWorkspace_tocorr_pf[k].var('u_'+sBin)


            mc_uPValZlike_pf = invertCdf(fromrdf["pfuP"+str(k+1)+"_uncorrected"][j], toCorr_cdf_pf, toCorr_xuP_pf, mc_cdf_pf, mc_xuP_pf)
            dt_uPValZlike_pf = invertCdf(fromrdf["pfuP"+str(k+1)+"_uncorrected"][j], toCorr_cdf_pf, toCorr_xuP_pf, data_cdf_pf, data_xuP_pf)

            tordf["pfuP"+str(k+1)+corrected+syst_postfix][j] += dt_uPValZlike_pf - mc_uPValZlike_pf

    rdf_tosave = ROOT.RDF.MakeNumpyDataFrame(tordf)
    rdf_tosave = calculateMET(rdf_tosave, corrected, mettype="", syst_postfix=syst_postfix)
    rdf_tosave = calculateMET(rdf_tosave, corrected, mettype="pf", syst_postfix=syst_postfix)

    if not os.path.exists(outfile) or args.overwrite:
        rdf_tosave.Snapshot("ntuple", outfile, [
            "uP1"+corrected+syst_postfix, "uP2"+corrected+syst_postfix, "met"+corrected+syst_postfix, "metphi"+corrected+syst_postfix,
            "pfuP1"+corrected+syst_postfix, "pfuP2"+corrected+syst_postfix, "pfmet"+corrected+syst_postfix, "pfmetphi"+corrected+syst_postfix
        ])
        # rdf_tosave.Snapshot("ntuple", outfile, list(set(original_cols + [
        #     "uP1"+corrected+syst_postfix, "uP2"+corrected+syst_postfix, "met"+corrected+syst_postfix, "metphi"+corrected+syst_postfix,
        #     "pfuP1"+corrected+syst_postfix, "pfuP2"+corrected+syst_postfix, "pfmet"+corrected+syst_postfix, "pfmetphi"+corrected+syst_postfix
        # ])))

    print(">>> "+infile+" success")


if __name__=='__main__':
    args = parse_args()

    ws_dicts = [
        {
            "data_mm":  "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_data_triple_muon_sigAndBck/", 
            "Zmm":      "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_WposMC_triple_muon_sigOnly/",
            "data_mm_pf":  "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_data_triple_muon_sigAndBck/", 
            "Zmm_pf":      "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg_pf":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos_pf":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_WposMC_triple_muon_sigOnly/",
        },

        {
            "data_mm":  "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_data_double_muon_sigAndBck/", 
            "Zmm":      "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_ZllMC_double_muon_sigOnly/",
            "Wmneg":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_WnegMC_double_muon_sigOnly/",
            "Wmpos":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_WposMC_double_muon_sigOnly/",
            "data_mm_pf":  "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_data_double_muon_sigAndBck/", 
            "Zmm_pf":      "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_ZllMC_double_muon_sigOnly/",
            "Wmneg_pf":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_WnegMC_double_muon_sigOnly/",
            "Wmpos_pf":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_WposMC_double_muon_sigOnly/",
        },

        {
            "data_mm":  "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_data_triple_muon_sigOnly/", 
            "Zmm":      "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/met_uncorrected_WposMC_triple_muon_sigOnly/",
            "data_mm_pf":  "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_data_triple_muon_sigOnly/", 
            "Zmm_pf":      "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg_pf":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos_pf":    "/home/moh/CROWN_working/Z_early_Run3/corrections/recoil/KitRecoilCorrections/Run3V02_outputs/pfmet_uncorrected_WposMC_triple_muon_sigOnly/",
        },
    ]

    ntuples = glob.glob(args.inpath)
    arguments = [{"ntuple": ntuple, "args": args, "ws_dict": ws_dict} for ntuple in ntuples for ws_dict in ws_dicts]

    # HERE tmp
    # tmp = []
    # for x in arguments:
    #     if "DYtoLL" in x["ntuple"]:
    #         tmp.append(x)
    #         break
    # arguments = tmp

    nthreads = 64
    if nthreads > len(arguments):
        nthreads = len(arguments)

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
