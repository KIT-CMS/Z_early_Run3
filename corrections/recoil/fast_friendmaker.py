#from sys import base_prefix
import ROOT
import os, sys
import gc
from tqdm import tqdm
import numpy as np
import glob
from multiprocessing import Pool, current_process, RLock

ROOT.gSystem.Load("./utils/PdfDiagonalizer_cc.so")

def load_fit_results(pathName, nBins, npars1=None, npars2=None):
    """
    function to get results from fit to mc and data z hadronic recoil

    (str) pathName: path to root file with RooWorkspace
    (int) nBins: 
    (int) npars1: 
    (int) npars2:
    """
    # set ROOT options to be silent
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.gROOT.ProcessLine(
        "RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);"
        )
    ROOT.gROOT.ProcessLine(
        "RooMsgService::instance().setSilentMode(true);"
        )

    # loop over both components of hadronic recoil and load fit results to ws
    ws = []

    for i in ["1", "2"]:
        file = ROOT.TFile(pathName+"pdfsU" + i + ".root")
        ws.append(file.Get("pdfsU"+i))

        # loop over number of Z pT bins and calculate cdfs
        for j in range(nBins):
            pdf = ws[-1].pdf("sig_" + str(j))
            myX = ws[-1].var("u_" + str(j))
            result = file.Get("fitResultU"+str(i)+"_"+str(j))
            cdf = pdf.createCdf(myX)
            ws[-1].Import(
                cdf, 
                ROOT.RooFit.RecycleConflictNodes(), 
                ROOT.RooFit.Silence()
                )

            # for stat unc
            npars = result.floatParsFinal().getSize()
            if i == "1":
                if npars1 == None:
                    npars1 = npars
                else:
                    assert npars1 == npars
            if i == "2":
                if npars2 == None:
                    npars2 = npars
                else:
                    assert npars2 == npars

            for ipar in range(npars):
                pdf_var_up = None
                pdf_var_dn = None
                cdf_var_up = None
                cdf_var_dn = None
                if result.status() != 0 or result.covQual() !=3:
                    pdf_var_up = pdf.Clone(pdf.GetName()+"_"+f"eig_{j}_ipar_{ipar}_Up")
                    pdf_var_dn = pdf.Clone(pdf.GetName()+"_"+f"eig_{j}_ipar_{ipar}_Down")
                    cdf_var_up = pdf_var_up.createCdf(myX)
                    cdf_var_dn = pdf_var_dn.createCdf(myX)
                else:
                    diago_up = ROOT.PdfDiagonalizer(f"eig_{j}_ipar_{ipar}_Up", ws[-1], result)
                    diago_dn = ROOT.PdfDiagonalizer(f"eig_{j}_ipar_{ipar}_Down", ws[-1], result)
                    pdf_dia_up = diago_up.diagonalize(pdf)
                    pdf_dia_dn = diago_dn.diagonalize(pdf)
                    pdf_var_up = diago_up.diagonalizeWithEigenVariations(pdf_dia_up, result, ipar, 1)
                    pdf_var_dn = diago_dn.diagonalizeWithEigenVariations(pdf_dia_dn, result, ipar, -1)

                    valid_up = True
                    valid_dn = True
                    for x in np.linspace(myX.getMin()+1, myX.getMax()-1, num=int(myX.getMax()-myX.getMin())):
                        myX.setVal(x)
                        val_up = pdf_var_up.getVal()
                        if val_up < 0. or val_up > 1. or np.isnan(val_up):
                            # pdf_var_up = diago_up.diagonalizeWithEigenVariations(pdf_dia_up, result, ipar, 0)
                            valid_up = False
                            break

                    for x in np.linspace(myX.getMin()+1, myX.getMax()-1, num=int(myX.getMax()-myX.getMin())):
                        myX.setVal(x)
                        val_dn = pdf_var_dn.getVal()
                        if val_dn < 0. or val_dn > 1. or np.isnan(val_dn):
                            # pdf_var_dn = diago_dn.diagonalizeWithEigenVariations(pdf_dia_dn, result, ipar, 0)
                            valid_dn = False
                            break

                    cdf_var_up = pdf_var_up.createCdf(myX)
                    for x in np.linspace(myX.getMin()+1, myX.getMax()-1, num=int(myX.getMax()-myX.getMin())):
                        myX.setVal(x)
                        val_up = cdf_var_up.getVal()
                        if val_up < 0. or val_up > 1. or np.isnan(val_up):
                            # pdf_var_up = diago_up.diagonalizeWithEigenVariations(pdf_dia_up, result, ipar, 0)
                            # cdf_var_up = pdf_var_up.createCdf(myX)
                            valid_up = False
                            break

                    cdf_var_dn = pdf_var_dn.createCdf(myX)
                    for x in np.linspace(myX.getMin()+1, myX.getMax()-1, num=int(myX.getMax()-myX.getMin())):
                        myX.setVal(x)
                        val_dn = cdf_var_dn.getVal()
                        if val_dn < 0. or val_dn > 1. or np.isnan(val_dn):
                            # pdf_var_dn = diago_dn.diagonalizeWithEigenVariations(pdf_dia_dn, result, ipar, 0)
                            # cdf_var_dn = pdf_var_dn.createCdf(myX)
                            valid_dn = False
                            break

                    if not (valid_up and valid_dn):
                        pdf_var_up = pdf.Clone(pdf.GetName()+"_"+f"eig_{j}_ipar_{ipar}_Up")
                        pdf_var_dn = pdf.Clone(pdf.GetName()+"_"+f"eig_{j}_ipar_{ipar}_Down")
                        cdf_var_up = pdf_var_up.createCdf(myX)
                        cdf_var_dn = pdf_var_dn.createCdf(myX)

                myX.setVal(0.5*(myX.getMax()+myX.getMin()))
                assert (cdf_var_up.getVal() > 0. and cdf_var_up.getVal() < 1.), f"{cdf_var_up.getVal()}, {cdf_var_up.getVal()}"
                assert (cdf_var_dn.getVal() > 0. and cdf_var_dn.getVal() < 1.), f"{cdf_var_dn.getVal()}, {cdf_var_dn.getVal()}"

                ws[-1].Import(cdf_var_up, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
                ws[-1].Import(cdf_var_dn, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
        file.Close()
    return ws, npars1, npars2


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

    xy = ROOT.TChain("ntuple")
    lep = ROOT.TChain("ntuple")
    xy.Add(inpath.replace("ntuples", "friends/xy"))
    lep.Add(inpath.replace("ntuples", "friends/lepton"))

    chain.AddFriend(xy)
    chain.AddFriend(lep)

    rdf = ROOT.RDataFrame(chain)

    return rdf, chain


def calculateMET(rdf, corrected, data=False, mettype="", postfix=""):
    rdf = rdf.Define(mettype+"metPx"+corrected+postfix, "-pt_vis_corr*cos(phi_vis_corr) - "+mettype+"uP1"+corrected+postfix+"*cos(bosonphi) + "+mettype+"uP2"+corrected+postfix+"*sin(bosonphi)")
    rdf = rdf.Define(mettype+"metPy"+corrected+postfix, "-pt_vis_corr*sin(phi_vis_corr) - "+mettype+"uP1"+corrected+postfix+"*sin(bosonphi) - "+mettype+"uP2"+corrected+postfix+"*cos(bosonphi)")

    rdf = rdf.Define(mettype+"met"+corrected+postfix, "sqrt("+mettype+"metPx"+corrected+postfix+"*"+mettype+"metPx"+corrected+postfix+" + "+mettype+"metPy"+corrected+postfix+"*"+mettype+"metPy"+corrected+postfix+")")
    rdf = rdf.Define(mettype+"metphi"+corrected+postfix, "atan2("+mettype+"metPy"+corrected+postfix+", "+mettype+"metPx"+corrected+postfix+")")

    if postfix=="":
        for updn in ["_up", "_dn"]:
            rdf = rdf.Define(
                mettype+"metPx"+corrected+updn, 
                "-pt_vis_corr"+updn+"*cos(phi_vis_corr"+updn+") - "+mettype+"uP1"+corrected+"*cos(bosonphi) + "+mettype+"uP2"+corrected+"*sin(bosonphi)"
                )
            rdf = rdf.Define(
                mettype+"metPy"+corrected+updn, 
                "-pt_vis_corr"+updn+"*sin(phi_vis_corr"+updn+") - "+mettype+"uP1"+corrected+"*sin(bosonphi) - "+mettype+"uP2"+corrected+"*cos(bosonphi)"
                )

            rdf = rdf.Define(
                mettype+"met"+corrected+updn, 
                "sqrt("+mettype+"metPx"+corrected+updn+"*"+mettype+"metPx"+corrected+updn+" + "+mettype+"metPy"+corrected+updn+"*"+mettype+"metPy"+corrected+updn+")"
                )
            rdf = rdf.Define(
                mettype+"metphi"+corrected+updn, 
                "atan2("+mettype+"metPy"+corrected+updn+", "+mettype+"metPx"+corrected+updn+")"
                )
    return rdf


def run(input_dict):
    infile = input_dict["ntuple"]
    ws_dict = input_dict["ws_dict"]

    process = infile.split('/')[-3]
    channel = infile.split('/')[-2]


    syst_postfix = ""
    if "_sigOnly" in ws_dict["data_mm"]:
        syst_postfix = "_sigOnly"
    elif "_double" in ws_dict["data_mm"]:
        syst_postfix = "_double"
    elif "_CB" in ws_dict["data_mm"]:
        syst_postfix = "_CB"

    if "_zrap0" in ws_dict["data_mm"]:
        syst_postfix = "_zrap0"
    elif "_zrap1" in ws_dict["data_mm"]:
        syst_postfix = "_zrap1"
    elif "_zrap2" in ws_dict["data_mm"]:
        syst_postfix = "_zrap2"

    doStatUnc = (syst_postfix == "")

    outfile = infile.replace("/ntuples/", "/friends/met"+syst_postfix+"/")
        

    outdir = outfile.replace(outfile.split('/')[-1], "")

    # HERE
    if os.path.isfile(outfile):
        f_tmp = None
        try:
            f_tmp = ROOT.TFile(outfile)
        except OSError:
            os.remove(outfile)
        if f_tmp and not f_tmp.IsZombie():
            return 1    
    # return -1
    # print(outdir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        # print("Created new directory: ", outdir)


    zPtBinEdges = [0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000]
    nBins = len(zPtBinEdges)-1

    corr = "_corr"

    rdf, _chain = prep(infile)


    if channel == 'mmet':
        # dataWorkspace, npars1, npars2 = load_fit_results(ws_dict["data_mm"], nBins)
        # mcWorkspace, npars1, npars2 = load_fit_results(ws_dict["Zmm"], nBins, npars1, npars2)
        # mcWorkspace_pos_tocorr, npars1, npars2 = load_fit_results(ws_dict["Wmpos"], nBins, npars1, npars2)
        # mcWorkspace_neg_tocorr, npars1, npars2 = load_fit_results(ws_dict["Wmneg"], nBins, npars1, npars2)
        dataWorkspace_pf, npars1, npars2 = load_fit_results(ws_dict["data_mm_pf"], nBins)
        mcWorkspace_pf, npars1, npars2 = load_fit_results(ws_dict["Zmm_pf"], nBins, npars1, npars2)
        mcWorkspace_pos_tocorr_pf, npars1, npars2 = load_fit_results(ws_dict["Wmpos_pf"], nBins, npars1, npars2)
        mcWorkspace_neg_tocorr_pf, npars1, npars2 = load_fit_results(ws_dict["Wmneg_pf"], nBins, npars1, npars2)
    elif channel == 'mm':
        # dataWorkspace, npars1, npars2 = load_fit_results(ws_dict["data_mm"], nBins)
        # mcWorkspace, npars1, npars2 = load_fit_results(ws_dict["Zmm"], nBins, npars1, npars2)
        # mcWorkspace_tocorr, npars1, npars2 = load_fit_results(ws_dict["Zmm"], nBins, npars1, npars2)
        dataWorkspace_pf, npars1, npars2 = load_fit_results(ws_dict["data_mm_pf"], nBins)
        mcWorkspace_pf, npars1, npars2 = load_fit_results(ws_dict["Zmm_pf"], nBins, npars1, npars2)
        mcWorkspace_tocorr_pf, npars1, npars2 = load_fit_results(ws_dict["Zmm_pf"], nBins, npars1, npars2)
    else:
        print('channel cannot be matched to mm, ee, mmet, emet')
        return -1
    

    # for other processes than DY: just use uncorrected values
    if not ("DYto2L" in process or "DYtoLL" in process or "WtoLNu" in process):
        # rdf = rdf.Define("uP1"+corr+syst_postfix, "uP1_uncorrected")
        # rdf = rdf.Define("uP2"+corr+syst_postfix, "uP2_uncorrected")

        rdf = rdf.Define("pfuP1"+corr+syst_postfix, "pfuP1_uncorrected")
        rdf = rdf.Define("pfuP2"+corr+syst_postfix, "pfuP2_uncorrected")

        # rdf = calculateMET(rdf, corr, mettype="", postfix=syst_postfix)
        rdf = calculateMET(rdf, corr, mettype="pf", postfix=syst_postfix)

        cols_stat_unc = []
        if doStatUnc:
            for stat_unc_idx in range(npars1+npars2):
                for updn in ["Up", "Down"]:
                    # rdf = rdf.Define("uP1"+corr+f"_stat{stat_unc_idx}{updn}", "uP1_uncorrected")
                    # rdf = rdf.Define("uP2"+corr+f"_stat{stat_unc_idx}{updn}", "uP2_uncorrected")

                    rdf = rdf.Define("pfuP1"+corr+f"_stat{stat_unc_idx}{updn}", "pfuP1_uncorrected")
                    rdf = rdf.Define("pfuP2"+corr+f"_stat{stat_unc_idx}{updn}", "pfuP2_uncorrected")
                
                    # rdf = calculateMET(rdf, corr, mettype="", postfix=f"_stat{stat_unc_idx}{updn}")
                    rdf = calculateMET(rdf, corr, mettype="pf", postfix=f"_stat{stat_unc_idx}{updn}")

                    # cols_stat_unc.append("met"+corr+f"_stat{stat_unc_idx}{updn}")
                    # cols_stat_unc.append("metphi"+corr+f"_stat{stat_unc_idx}{updn}")
                    cols_stat_unc.append("pfmet"+corr+f"_stat{stat_unc_idx}{updn}")
                    cols_stat_unc.append("pfmetphi"+corr+f"_stat{stat_unc_idx}{updn}")
        
        cols_lepunc = []
        if syst_postfix=="":
            cols_lepunc += [
                #"met_corr_up", "met_corr_dn", "metphi_corr_up", "metphi_corr_dn",
                "pfmet_corr_up", "pfmet_corr_dn", "pfmetphi_corr_up", "pfmetphi_corr_dn"
            ]

        rdf.Snapshot("ntuple", outfile, [
            #"uP1"+corr+syst_postfix, "uP2"+corr+syst_postfix, 
            #"met"+corr+syst_postfix, "metphi"+corr+syst_postfix,
            "pfuP1"+corr+syst_postfix, "pfuP2"+corr+syst_postfix, 
            "pfmet"+corr+syst_postfix, "pfmetphi"+corr+syst_postfix
        ] + cols_stat_unc + cols_lepunc)
            # rdf.Snapshot("ntuple", outfile, list(set(original_cols + [
            #     "uP1"+corrected+syst_postfix, "uP2"+corrected+syst_postfix, "met"+corrected+syst_postfix, "metphi"+corrected+syst_postfix,
            #     "pfuP1"+corrected+syst_postfix, "pfuP2"+corrected+syst_postfix, "pfmet"+corrected+syst_postfix, "pfmetphi"+corrected+syst_postfix
            # ])))

        if channel == 'mmet':
            #del dataWorkspace
            #del mcWorkspace
            #del mcWorkspace_pos_tocorr
            #del mcWorkspace_neg_tocorr
            del dataWorkspace_pf
            del mcWorkspace_pf
            del mcWorkspace_pos_tocorr_pf
            del mcWorkspace_neg_tocorr_pf
        elif channel == 'mm':
            #del dataWorkspace
            #del mcWorkspace
            #del mcWorkspace_tocorr
            del dataWorkspace_pf
            del mcWorkspace_pf
            del mcWorkspace_tocorr_pf
        del rdf
        gc.collect()

        return 1





    # Actual correction
    # rdf = rdf.Define("uP1"+corr+syst_postfix, "uP1_uncorrected")
    # rdf = rdf.Define("uP2"+corr+syst_postfix, "uP2_uncorrected")
    rdf = rdf.Define("pfuP1"+corr+syst_postfix, "pfuP1_uncorrected")
    rdf = rdf.Define("pfuP2"+corr+syst_postfix, "pfuP2_uncorrected")
    cols_stat_unc_u = []
    if doStatUnc:
        for updn in ["Up", "Down"]:
            for stat_unc_idx in range(npars1+npars2):
                # rdf = rdf.Define("uP1"+corr+f"_stat{stat_unc_idx}{updn}", "uP1_uncorrected")
                # rdf = rdf.Define("uP2"+corr+f"_stat{stat_unc_idx}{updn}", "uP2_uncorrected")
                rdf = rdf.Define("pfuP1"+corr+f"_stat{stat_unc_idx}{updn}", "pfuP1_uncorrected")
                rdf = rdf.Define("pfuP2"+corr+f"_stat{stat_unc_idx}{updn}", "pfuP2_uncorrected")
                # cols_stat_unc_u.append("uP1"+corr+f"_stat{stat_unc_idx}{updn}")
                # cols_stat_unc_u.append("uP2"+corr+f"_stat{stat_unc_idx}{updn}")
                cols_stat_unc_u.append("pfuP1"+corr+f"_stat{stat_unc_idx}{updn}")
                cols_stat_unc_u.append("pfuP2"+corr+f"_stat{stat_unc_idx}{updn}")

    columns = [
        # "uP1"+corr+syst_postfix,
        # "uP2"+corr+syst_postfix,
        # "uP1_uncorrected",
        # "uP2_uncorrected",
        "pfuP1"+corr+syst_postfix,
        "pfuP2"+corr+syst_postfix,
        "pfuP1_uncorrected",
        "pfuP2_uncorrected",
        "pt_vis_corr", "pt_vis_corr_up", "pt_vis_corr_dn",
        "phi_vis_corr", "phi_vis_corr_up", "phi_vis_corr_dn",
        "bosonpt",
        "bosonphi",
        'q_1'
    ] + cols_stat_unc_u
    fromrdf = rdf.AsNumpy(columns = columns)
    tordf = {
        # "uP1"+corr+syst_postfix: fromrdf["uP1"+corr+syst_postfix],
        # "uP2"+corr+syst_postfix: fromrdf["uP2"+corr+syst_postfix],
        # "uP1_uncorrected": fromrdf["uP1_uncorrected"],
        # "uP2_uncorrected": fromrdf["uP2_uncorrected"],
        "pfuP1"+corr+syst_postfix: fromrdf["pfuP1"+corr+syst_postfix],
        "pfuP2"+corr+syst_postfix: fromrdf["pfuP2"+corr+syst_postfix],
        "pfuP1_uncorrected": fromrdf["pfuP1_uncorrected"],
        "pfuP2_uncorrected": fromrdf["pfuP2_uncorrected"],
        "pt_vis_corr": fromrdf["pt_vis_corr"],
        "pt_vis_corr_up": fromrdf["pt_vis_corr_up"],
        "pt_vis_corr_dn": fromrdf["pt_vis_corr_dn"],
        "phi_vis_corr": fromrdf["phi_vis_corr"],
        "phi_vis_corr_up": fromrdf["phi_vis_corr_up"],
        "phi_vis_corr_dn": fromrdf["phi_vis_corr_dn"],
        "bosonpt": fromrdf["bosonpt"],
        "bosonphi": fromrdf["bosonphi"]
    }
    if doStatUnc:
        for col in cols_stat_unc_u:
            tordf[col] = fromrdf[col]

    for j in range(len(fromrdf["bosonpt"])):
        # Get pT bin of genboson
        kBin = nBins-1
        for k in range(nBins):
            if ((channel == 'mm' or channel == 'ee') and fromrdf["pt_vis_corr"][j] < zPtBinEdges[k]) or ((channel == 'mmet' or channel == 'emet') and fromrdf["bosonpt"][j] < zPtBinEdges[k]): # TODO check if right order
                kBin = k-1
                break
        sBin = str(kBin)

        # MET corrections for parallel and perpendicular recoil component
        for k in [0,1]:
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

            tordf["uP"+str(k+1)+corr+syst_postfix][j] += dt_uPValZlike - mc_uPValZlike
            """

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

            tordf["pfuP"+str(k+1)+corr+syst_postfix][j] += dt_uPValZlike_pf - mc_uPValZlike_pf

            if doStatUnc:
                for updn in ["Up", "Down"]:
                    npars = npars1 if k == 0 else npars2
                    stat_unc_idx = 0
                    for ipar in range(npars):
                        """
                        # PUPPI
                        data_cdf   = dataWorkspace[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                        data_xuP   = dataWorkspace[k].var('u_'+sBin)

                        if ('met' in channel) and (fromrdf["q_1"][j] == 1):
                            mc_cdf     = mcWorkspace[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            toCorr_cdf = mcWorkspace_pos_tocorr[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            mc_xuP     = mcWorkspace[k].var('u_'+sBin)
                            toCorr_xuP = mcWorkspace_pos_tocorr[k].var('u_'+sBin)
                        elif ('met' in channel) and (fromrdf["q_1"][j] == -1):
                            mc_cdf     = mcWorkspace[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            toCorr_cdf = mcWorkspace_neg_tocorr[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            mc_xuP     = mcWorkspace[k].var('u_'+sBin)
                            toCorr_xuP = mcWorkspace_neg_tocorr[k].var('u_'+sBin)
                        else:
                            mc_cdf     = mcWorkspace[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            toCorr_cdf = mcWorkspace_tocorr[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            mc_xuP     = mcWorkspace[k].var('u_'+sBin)
                            toCorr_xuP = mcWorkspace_tocorr[k].var('u_'+sBin)

                        mc_uPValZlike = invertCdf(fromrdf["uP"+str(k+1)+"_uncorrected"][j], toCorr_cdf, toCorr_xuP, mc_cdf, mc_xuP)
                        dt_uPValZlike = invertCdf(fromrdf["uP"+str(k+1)+"_uncorrected"][j], toCorr_cdf, toCorr_xuP, data_cdf, data_xuP)

                        tordf["uP"+str(k+1)+corr+f"_stat{stat_unc_idx}{updn}"][j] += dt_uPValZlike - mc_uPValZlike
                        """
                        # PF
                        data_cdf_pf   = dataWorkspace_pf[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                        data_xuP_pf   = dataWorkspace_pf[k].var('u_'+sBin)

                        if ('met' in channel) and (fromrdf["q_1"][j] == 1):
                            mc_cdf_pf     = mcWorkspace_pf[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            toCorr_cdf_pf = mcWorkspace_pos_tocorr_pf[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            mc_xuP_pf     = mcWorkspace_pf[k].var('u_'+sBin)
                            toCorr_xuP_pf = mcWorkspace_pos_tocorr_pf[k].var('u_'+sBin)
                        elif ('met' in channel) and (fromrdf["q_1"][j] == -1):
                            mc_cdf_pf     = mcWorkspace_pf[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            toCorr_cdf_pf = mcWorkspace_neg_tocorr_pf[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            mc_xuP_pf     = mcWorkspace_pf[k].var('u_'+sBin)
                            toCorr_xuP_pf = mcWorkspace_neg_tocorr_pf[k].var('u_'+sBin)
                        else:
                            mc_cdf_pf     = mcWorkspace_pf[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            toCorr_cdf_pf = mcWorkspace_tocorr_pf[k].function(f"sig_{sBin}_eig_{sBin}_ipar_{ipar}_{updn}_cdf_Int[u_{sBin}_prime|CDF]_Norm[u_{sBin}_prime]")
                            mc_xuP_pf     = mcWorkspace_pf[k].var('u_'+sBin)
                            toCorr_xuP_pf = mcWorkspace_tocorr_pf[k].var('u_'+sBin)

                        mc_uPValZlike_pf = invertCdf(fromrdf["pfuP"+str(k+1)+"_uncorrected"][j], toCorr_cdf_pf, toCorr_xuP_pf, mc_cdf_pf, mc_xuP_pf)
                        dt_uPValZlike_pf = invertCdf(fromrdf["pfuP"+str(k+1)+"_uncorrected"][j], toCorr_cdf_pf, toCorr_xuP_pf, data_cdf_pf, data_xuP_pf)

                        tordf["pfuP"+str(k+1)+corr+f"_stat{stat_unc_idx}{updn}"][j] += dt_uPValZlike_pf - mc_uPValZlike_pf

                        stat_unc_idx += 1

    rdf_tosave = ROOT.RDF.MakeNumpyDataFrame(tordf)
    # rdf_tosave = calculateMET(rdf_tosave, corr, mettype="", postfix=syst_postfix)
    rdf_tosave = calculateMET(rdf_tosave, corr, mettype="pf", postfix=syst_postfix)

    cols_stat_unc = []
    if doStatUnc:
        for updn in ["Up", "Down"]:
            for stat_unc_idx in range(npars1+npars2):
                # rdf_tosave = calculateMET(rdf_tosave, corr, mettype="", postfix=f"_stat{stat_unc_idx}{updn}")
                rdf_tosave = calculateMET(rdf_tosave, corr, mettype="pf", postfix=f"_stat{stat_unc_idx}{updn}")
                # cols_stat_unc.append("met"+corr+f"_stat{stat_unc_idx}{updn}")
                # cols_stat_unc.append("metphi"+corr+f"_stat{stat_unc_idx}{updn}")
                cols_stat_unc.append("pfmet"+corr+f"_stat{stat_unc_idx}{updn}")
                cols_stat_unc.append("pfmetphi"+corr+f"_stat{stat_unc_idx}{updn}")
    
    cols_lepunc = []
    if syst_postfix=="":
        cols_lepunc += [
            # "met_corr_up", "met_corr_dn", "metphi_corr_up", "metphi_corr_dn",
            "pfmet_corr_up", "pfmet_corr_dn", "pfmetphi_corr_up", "pfmetphi_corr_dn"
        ]

    rdf_tosave.Snapshot("ntuple", outfile, [
        # "uP1"+corr+syst_postfix, "uP2"+corr+syst_postfix, "met"+corr+syst_postfix, "metphi"+corr+syst_postfix,
        "pfuP1"+corr+syst_postfix, "pfuP2"+corr+syst_postfix, "pfmet"+corr+syst_postfix, "pfmetphi"+corr+syst_postfix
    ] + cols_stat_unc + cols_lepunc)

    if channel == 'mmet':
        #del dataWorkspace
        #del mcWorkspace
        #del mcWorkspace_pos_tocorr
        #del mcWorkspace_neg_tocorr
        del dataWorkspace_pf
        del mcWorkspace_pf
        del mcWorkspace_pos_tocorr_pf
        del mcWorkspace_neg_tocorr_pf
    elif channel == 'mm':
        #del dataWorkspace
        #del mcWorkspace
        #del mcWorkspace_tocorr
        del dataWorkspace_pf
        del mcWorkspace_pf
        del mcWorkspace_tocorr_pf
    del fromrdf
    del rdf_tosave
    del rdf
    gc.collect()

    return 1


if __name__=='__main__':
    ROOT.gROOT.SetBatch(True) 

    infiles = '/ceph/jdriesch/CROWN_samples/Run3V07/ntuples/2022/*/mm*/*.root'

    basepath = '/work/jdriesch/earlyrun3/Z_early_Run3/corrections/recoil/'\
        'KitRecoilCorrections/Run3V07_ptuncorr_outputs/'

    ws_dicts = [
        {
            "data_mm": basepath+"met_uncorrected_data_triple_muon_sigAndBck/", 
            "Zmm": basepath+"met_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg": basepath+"met_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos": basepath+"met_uncorrected_WposMC_triple_muon_sigOnly/",
            "data_mm_pf": basepath+"pfmet_uncorrected_data_triple_muon_sigAndBck/", 
            "Zmm_pf": basepath+"pfmet_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg_pf": basepath+"pfmet_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos_pf": basepath+"pfmet_uncorrected_WposMC_triple_muon_sigOnly/",
        },

        {
            "data_mm": basepath+"met_uncorrected_data_double_muon_sigAndBck/", 
            "Zmm": basepath+"met_uncorrected_ZllMC_double_muon_sigOnly/",
            "Wmneg": basepath+"met_uncorrected_WnegMC_double_muon_sigOnly/",
            "Wmpos": basepath+"met_uncorrected_WposMC_double_muon_sigOnly/",
            "data_mm_pf": basepath+"pfmet_uncorrected_data_double_muon_sigAndBck/", 
            "Zmm_pf": basepath+"pfmet_uncorrected_ZllMC_double_muon_sigOnly/",
            "Wmneg_pf": basepath+"pfmet_uncorrected_WnegMC_double_muon_sigOnly/",
            "Wmpos_pf": basepath+"pfmet_uncorrected_WposMC_double_muon_sigOnly/",
        },

        {
            "data_mm": basepath+"met_uncorrected_data_triple_muon_sigOnly/", 
            "Zmm": basepath+"met_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg": basepath+"met_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos": basepath+"met_uncorrected_WposMC_triple_muon_sigOnly/",
            "data_mm_pf": basepath+"pfmet_uncorrected_data_triple_muon_sigOnly/", 
            "Zmm_pf": basepath+"pfmet_uncorrected_ZllMC_triple_muon_sigOnly/",
            "Wmneg_pf": basepath+"pfmet_uncorrected_WnegMC_triple_muon_sigOnly/",
            "Wmpos_pf": basepath+"pfmet_uncorrected_WposMC_triple_muon_sigOnly/",
        },
         
        {
            "data_mm": basepath+"met_uncorrected_data_triple_muon_sigAndBck_zrap0/",
            "Zmm": basepath+"met_uncorrected_ZllMC_triple_muon_sigOnly_zrap0/",
            "Wmneg": basepath+"met_uncorrected_WnegMC_triple_muon_sigOnly_zrap0/",
            "Wmpos": basepath+"met_uncorrected_WposMC_triple_muon_sigOnly_zrap0/",
            "data_mm_pf": basepath+"pfmet_uncorrected_data_triple_muon_sigAndBck_zrap0/",
            "Zmm_pf": basepath+"pfmet_uncorrected_ZllMC_triple_muon_sigOnly_zrap0/",
            "Wmneg_pf": basepath+"pfmet_uncorrected_WnegMC_triple_muon_sigOnly_zrap0/",
            "Wmpos_pf": basepath+"pfmet_uncorrected_WposMC_triple_muon_sigOnly_zrap0/",
        },

        {
            "data_mm": basepath+"met_uncorrected_data_triple_muon_sigAndBck_zrap1/",
            "Zmm": basepath+"met_uncorrected_ZllMC_triple_muon_sigOnly_zrap1/",
            "Wmneg": basepath+"met_uncorrected_WnegMC_triple_muon_sigOnly_zrap1/",
            "Wmpos": basepath+"met_uncorrected_WposMC_triple_muon_sigOnly_zrap1/",
            "data_mm_pf": basepath+"pfmet_uncorrected_data_triple_muon_sigAndBck_zrap1/",
            "Zmm_pf": basepath+"pfmet_uncorrected_ZllMC_triple_muon_sigOnly_zrap1/",
            "Wmneg_pf": basepath+"pfmet_uncorrected_WnegMC_triple_muon_sigOnly_zrap1/",
            "Wmpos_pf": basepath+"pfmet_uncorrected_WposMC_triple_muon_sigOnly_zrap1/",
        },

        {
            "data_mm": basepath+"met_uncorrected_data_triple_muon_sigAndBck_zrap2/",
            "Zmm": basepath+"met_uncorrected_ZllMC_triple_muon_sigOnly_zrap2/",
            "Wmneg": basepath+"met_uncorrected_WnegMC_triple_muon_sigOnly_zrap2/",
            "Wmpos": basepath+"met_uncorrected_WposMC_triple_muon_sigOnly_zrap2/",
            "data_mm_pf": basepath+"pfmet_uncorrected_data_triple_muon_sigAndBck_zrap2/",
            "Zmm_pf": basepath+"pfmet_uncorrected_ZllMC_triple_muon_sigOnly_zrap2/",
            "Wmneg_pf": basepath+"pfmet_uncorrected_WnegMC_triple_muon_sigOnly_zrap2/",
            "Wmpos_pf": basepath+"pfmet_uncorrected_WposMC_triple_muon_sigOnly_zrap2/",
        }, 
    ]

    ntuples = glob.glob(infiles)
    arguments = [{"ntuple": ntuple, "ws_dict": ws_dict} for ntuple in ntuples for ws_dict in ws_dicts]

    nthreads = 32
    with Pool(nthreads, initargs=(RLock(),), initializer=tqdm.set_lock, maxtasksperchild=1) as pool:
        for _ in tqdm(
            pool.imap_unordered(run, arguments, chunksize=1),
            total=len(arguments),
            desc="Total progess",
            # position=nthreads + 1,
            dynamic_ncols=True,
            leave=True,
        ):
            pass
        pool.close()
