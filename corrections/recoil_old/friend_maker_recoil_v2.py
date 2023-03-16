#from sys import base_prefix
import ROOT
import argparse
import os
from tqdm import tqdm
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="produce friend of input ntuple with corrected met")   
    #parser.add_argument('-W', '--wsdir', default='/work/jheitkoetter/MitEwk13TeV_CMSSW_94X/Recoil/completed_plots/plots_2018_fine_pf/', help='directory of RooFitWorkSpace')
    #parser.add_argument('-C', '--corr', default='Zmm', help='Correction Samples: Zmm, Wmpos, Wmneg, ...')
    parser.add_argument('-F', '--finalstate', default='mm', help='Final state: mm, ee, mmet or emet')
    parser.add_argument('-I', '--inpath', default='/ceph/moh/CROWN_samples/EarlyRun3_V08/ntuples/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/', help='path to input samples')
    parser.add_argument('--overwrite', action='store_true', default=False)
    parser.add_argument('--test', action='store_true', default=False)   
    args = parser.parse_args()
    return args

def load_fit_results(pathName, nBins):
    # loop over parallel and perpendicular component of hadronic recoil and load fit results

    workspace = []

    for i in ["1","2"]:
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
    rdf = ROOT.RDataFrame("ntuple", inpath)
    
    # check how many leptons in final state and adjust total leptonic momentum accordingly
    finalstate = inpath.split("/")[-2]

    # apply rochester corrections if muon involved (mm or mmet)
    if "mm" in finalstate:
        cor = "_rc"
    else:
        cor = ""

    if not ("DY" in inpath or "WJet" in inpath):
        bosonphi = "phi_vis_c"
        bosonpt = "pt_vis_c"
    else:
        bosonphi = "genbosonphi"
        bosonpt = "genbosonpt"

    if finalstate == "ee" or finalstate == "mm":
        rdf = rdf.Define("pt_vis_c_x", "pt"+cor+"_1*cos(phi_1) + pt"+cor+"_2*cos(phi_2)")
        rdf = rdf.Define("pt_vis_c_y", "pt"+cor+"_1*sin(phi_1) + pt"+cor+"_2*sin(phi_2)")
        rdf = rdf.Define("pt_vis_c", "sqrt(pt_vis_c_x*pt_vis_c_x + pt_vis_c_y*pt_vis_c_y)")
        rdf = rdf.Define("phi_vis_c", "atan2(pt_vis_c_y, pt_vis_c_x)")

    else:
        rdf = rdf.Define("pt_vis_c", "pt"+cor+"_1")
        rdf = rdf.Define("phi_vis_c", "phi_1")

    rdf = rdf.Define("uPx", met+"_uncorrected*cos("+met+"phi_uncorrected) + pt_vis_c*cos(phi_vis_c)") # sign is wrong but consistently -> met is correct
    rdf = rdf.Define("uPy", met+"_uncorrected*sin("+met+"phi_uncorrected) + pt_vis_c*sin(phi_vis_c)")

    rdf = rdf.Define("uP1_uncorrected", "- (uPx*cos("+bosonphi+") + uPy*sin("+bosonphi+"))") # maybe use reco-phi instead
    rdf = rdf.Define("uP2_uncorrected", "uPx*sin("+bosonphi+") - uPy*cos("+bosonphi+")")
    
    return rdf


def calculateMET(rdf, corrected, data=False, met="met"):
    rdf = rdf.Define(met+"Px"+corrected, "-pt_vis_c*cos(phi_vis_c) - uP1"+corrected+"*cos(genbosonphi) + uP2"+corrected+"*sin(genbosonphi)") # TODO check if name already used
    rdf = rdf.Define(met+"Py"+corrected, "-pt_vis_c*sin(phi_vis_c) - uP1"+corrected+"*sin(genbosonphi) - uP2"+corrected+"*cos(genbosonphi)")

    rdf = rdf.Define(met+corrected, "sqrt("+met+"Px"+corrected+"*"+met+"Px"+corrected+" + "+met+"Py"+corrected+"*"+met+"Py"+corrected+")")
    rdf = rdf.Define(met+"phi"+corrected, "atan2("+met+"Py"+corrected+", "+met+"Px"+corrected+")")

    return rdf


def main():
    args = parse_args()

    ws_dict = {
        "data_mm":  "/work/jheitkoetter/MitEwk13TeV_CMSSW_94X/Recoil/V12_outputs/met_uncorrected_data_triple_muon_sigOnly/", 
        "data_ee":  "/work/jheitkoetter/MitEwk13TeV_CMSSW_94X/Recoil/V12_outputs/met_uncorrected_data_triple_elec_sigOnly/", 
        "Zmm":      "/work/jheitkoetter/MitEwk13TeV_CMSSW_94X/Recoil/V12_outputs/met_uncorrected_ZllMC_triple_muon_sigOnly/",
        "Zee":      "/work/jheitkoetter/MitEwk13TeV_CMSSW_94X/Recoil/V12_outputs/met_uncorrected_ZllMC_triple_elec_sigOnly/",
        "Wmneg":    "/work/jdriesch/phd/Z_early_Run3/corrections/recoil/KitRecoilCorrections/V12_outputs/met_uncorrected_WnegMC_triple_muon_sigOnly/",
        "Wmpos":    "/work/jheitkoetter/MitEwk13TeV_CMSSW_94X/Recoil/V12_outputs/met_uncorrected_WposMC_triple_muon_sigOnly/",
        "Wepos":    "/work/jdriesch/phd/Z_early_Run3/corrections/recoil/KitRecoilCorrections/V12_outputs/met_uncorrected_WposMC_triple_elec_sigOnly/",
        "Weneg":    "/work/jheitkoetter/MitEwk13TeV_CMSSW_94X/Recoil/V12_outputs/met_uncorrected_WnegMC_triple_elec_sigOnly/",
        }

    zPtBinEdges = [0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000]
    nBins = len(zPtBinEdges)-1


    if args.finalstate == 'mmet':
        dataWorkspace = load_fit_results(ws_dict["data_mm"], nBins)
        mcWorkspace_pos = load_fit_results(ws_dict["Wmpos"], nBins)
        mcWorkspace_neg = load_fit_results(ws_dict["Wmneg"], nBins)           
        mcWorkspace_pos_tocorr = load_fit_results(ws_dict["Wmpos"], nBins)
        mcWorkspace_neg_tocorr = load_fit_results(ws_dict["Wmneg"], nBins)   

    elif args.finalstate == 'emet':
        dataWorkspace = load_fit_results(ws_dict["data_ee"], nBins)
        mcWorkspace_pos = load_fit_results(ws_dict["Wepos"], nBins)
        mcWorkspace_neg = load_fit_results(ws_dict["Weneg"], nBins)        
        mcWorkspace_pos_tocorr = load_fit_results(ws_dict["Wepos"], nBins)
        mcWorkspace_neg_tocorr = load_fit_results(ws_dict["Weneg"], nBins)

    elif args.finalstate == 'mm':
        dataWorkspace = load_fit_results(ws_dict["data_mm"], nBins)
        mcWorkspace = load_fit_results(ws_dict["Zmm"], nBins)
        mcWorkspace_tocorr = load_fit_results(ws_dict["Zmm"], nBins)

    elif args.finalstate == 'ee':
        dataWorkspace = load_fit_results(ws_dict["data_ee"], nBins)
        mcWorkspace = load_fit_results(ws_dict["Zee"], nBins)        
        mcWorkspace_tocorr = load_fit_results(ws_dict["Zee"], nBins)

    else:
        print('finalstate cannot be matched to mm, ee, mmet, emet')
        return

    version = args.inpath.split('/')[4]
    process = args.inpath.split('/')[-2]
    year = args.inpath.split('/')[-3]

    inpath = args.inpath + "/" + args.finalstate + "/" + process

    
    if args.test:
        outdir = '/ceph/jdriesch/CROWN_samples/test/'
    else:
        outdir = '/ceph/jdriesch/CROWN_samples/' + version + '/friends/metcorr/' + year + '/' + process + "/" + args.finalstate + "/"

    # make sure that variables are named consistently (for Run2 met is already corrected)
    if "18" in year:
        corrected = "_corrected"
        print("run2 data")
    else:
        corrected = ""
        print("run3 data")
    met = 'met'
    """
    if "_pf" in args.wsdir:
        met = "pfmet"
    else: 
        met = "met"
    """


    i = 0    

    while os.path.exists(inpath + "_"+str(i)+".root"):
        infile = inpath + "_"+str(i)+".root"
        outfile= outdir + infile.split("/")[-1]

        if not os.path.exists(outdir):
            os.makedirs(outdir)
            print("Created new directory: ", outdir)
            
        rdf = prep(infile, met)

        # for other processes than DY: just use uncorrected values
        if not ("DY" in process or "WJet" in process):
            if i==0:
                print("No corrections will be applied")
            rdf = rdf.Define("uP1"+corrected, "uP1_uncorrected")
            rdf = rdf.Define("uP2"+corrected, "uP2_uncorrected")
            rdf = rdf.Define(met+corrected, 'met_uncorrected')
            rdf = rdf.Define(met+"phi"+corrected, 'metphi_uncorrected')

            if not os.path.exists(outfile) or args.overwrite:
                rdf.Snapshot("ntuple", outfile, ["uP1_uncorrected", "uP2_uncorrected", "uP1"+corrected, "uP2"+corrected, met+corrected, met+"phi"+corrected, "pt_vis_c", "phi_vis_c"])

            if i%10==0: 
                print("done {} of {}".format(i, 100))

            i+=1

            continue


        rdf = rdf.Define("uP1"+corrected, "uP1_uncorrected")
        rdf = rdf.Define("uP2"+corrected, "uP2_uncorrected")


        columns = ["uP1_corrected", "uP2_corrected", "uP1_uncorrected", "uP2_uncorrected", "pt_vis_c", "phi_vis_c", "genbosonpt", "genbosonphi", 'q_1']
        fromrdf = rdf.AsNumpy(columns = columns)
        tordf = {
            "uP1"+corrected: fromrdf["uP1_corrected"], 
            "uP2"+corrected: fromrdf["uP2_corrected"],             
            "uP1_uncorrected": fromrdf["uP1_uncorrected"], 
            "uP2_uncorrected": fromrdf["uP2_uncorrected"], 
            "pt_vis_c": fromrdf["pt_vis_c"], 
            "phi_vis_c": fromrdf["phi_vis_c"], 
            "genbosonpt": fromrdf["genbosonpt"], 
            "genbosonphi": fromrdf["genbosonphi"]
            }

        for j in tqdm(range(len(fromrdf["genbosonpt"]))):
            # Get pT bin of genboson
            kBin = nBins-1
            for k in range(nBins):
                if ((args.finalstate == 'mm' or args.finalstate == 'ee') and fromrdf["pt_vis_c"][j] < zPtBinEdges[k]) or ((args.finalstate == 'mmet' or args.finalstate == 'emet') and fromrdf["genbosonpt"][j] < zPtBinEdges[k]): # TODO check if right order
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

                if ('met' in args.finalstate) and (fromrdf["q_1"][j] == 1):
                    mc_cdf     = mcWorkspace_pos[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')
                    toCorr_cdf = mcWorkspace_pos_tocorr[k].function('sig_'+sBin+'_cdf_Int[u_'+sBin+'_prime|CDF]_Norm[u_'+sBin+'_prime]')    
                    mc_xuP     = mcWorkspace_pos[k].var('u_'+sBin)
                    toCorr_xuP = mcWorkspace_pos_tocorr[k].var('u_'+sBin)

                elif ('met' in args.finalstate) and (fromrdf["q_1"][j] == -1):                   
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
                data_uPValZlike = invertCdf(mc_uPValZlike, mc_cdf, mc_xuP, data_cdf, data_xuP)

                tordf["uP"+str(k+1)+corrected][j] += data_uPValZlike - mc_uPValZlike

        rdf_tosave = ROOT.RDF.MakeNumpyDataFrame(tordf)
        rdf_tosave = calculateMET(rdf_tosave, corrected, met=met)     

        if not os.path.exists(outfile) or args.overwrite:
            rdf_tosave.Snapshot("ntuple", outfile, ["uP1_uncorrected", "uP2_uncorrected", "uP1"+corrected, "uP2"+corrected, met+corrected, met+"phi"+corrected, "pt_vis_c", "phi_vis_c"])

        if i%10==0: 
            print("done {} of {}".format(i, 100))

        i+=1


    print("great success")


if __name__=='__main__':
    main()