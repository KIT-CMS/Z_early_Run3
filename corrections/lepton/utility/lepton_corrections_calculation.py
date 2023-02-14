from ast import arg
import ROOT
import argparse
import os
from tqdm import tqdm
import numpy as np
import utils
import array as a
import json
import glob

def parse_args():
    parser = argparse.ArgumentParser(description="Lepton correction program. To obtain corrections, please proceed as follows: \n 1. Produce histograms by setting option -H\n 2. Calculate Z mass peak position and width via -C\n 3. Produce friend tree with corrected lepton pt values via -E\n 4. Plot corrected distributions via -P")
    parser.add_argument('-F', '--finalstate', default='mm', help='Final state: mm, ee, mmet or emet')
    parser.add_argument('-I', '--inpath', default='/ceph/moh/CROWN_samples/EarlyRun3_V08/ntuples/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/', help='path to input samples')
    parser.add_argument('--overwrite', action='store_true', default=False, help='Set to overwrite existing files')
    parser.add_argument('--test', action='store_true', default=False, help='set to put into test directory')
    parser.add_argument('-H', '--histogram', action = 'store_true', default=False, help='Create binned histograms')
    parser.add_argument('-C', '--calculation', action = 'store_true', default=False, help='Set histogram mode')
    parser.add_argument('-E', '--evaluation', action = 'store_true', default=False, help='Set histogram mode')
    parser.add_argument('-P', '--plot', action = 'store_true', default=False)
    parser.add_argument('--corr', action= 'store_true', default=False)
    parser.add_argument('--info', action= 'store_true', default=False)
    parser.add_argument('-R', '--run', default='2', help='Run number. Either 2 (condsiders 2018 only) or 3 (for 2022)')
    parser.add_argument('-V', '--version', default='v1')
    parser.add_argument('--finalize', action='store_true', default=False)
    parser.add_argument('-x', '--res-factor', default=1.)
    parser.add_argument('--run_i', default=355862)
    parser.add_argument('--run_f', default=357482) #357900
    args = parser.parse_args()
    return args


def get_paths(args, mode, kw = None):
    # get ntuple or friends paths
    if mode == 0:
        # paths with uncorrected ntuples
        if not args.corr:
            if args.run=='2':
                if args.finalstate=='mm':
                    mc = ['/ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/mm/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_*.root']
                    dt = ['/ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples/2018/SingleMuon_Run2018{run}-UL2018/mm/SingleMuon_Run2018{run}-UL2018_*.root'.format(run=r) for r in ['A', 'B', 'C', 'D']]         
                elif args.finalstate=='ee':
                    mc = ['/ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/ee/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_*.root']
                    dt = ['/ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples/2018/EGamma_Run2018{run}-UL2018/ee/EGamma_Run2018{run}-UL2018_*.root'.format(run=r) for r in ['A', 'B', 'C']]   
                    dt += ['/ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples/2018/EGamma_Run2018{run}-UL2018/ee/EGamma_Run2018{run}-UL2018_{i}.root'.format(run='D', i=i) for i in range(113)]
                    dt += ['/ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples/2018/EGamma_Run2018{run}-UL2018/ee/EGamma_Run2018{run}-UL2018_{i}.root'.format(run='D', i='113_ee')]
                    dt += ['/ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples/2018/EGamma_Run2018{run}-UL2018/ee/EGamma_Run2018{run}-UL2018_{i}.root'.format(run='D', i=i) for i in range(114, 117)]
                    dt += ['/ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples/2018/EGamma_Run2018{run}-UL2018/ee/EGamma_Run2018{run}-UL2018_{i}.root'.format(run='D', i='117_ee')]
                    dt += ['/ceph/moh/CROWN_samples/EarlyRun3_V12/ntuples/2018/EGamma_Run2018{run}-UL2018/ee/EGamma_Run2018{run}-UL2018_{i}.root'.format(run='D', i=i) for i in range(118, 173)] 
                else:
                    print('Only Files with two leptons of the same kind in the finalstate should be considered for corrections. Aborting...')
                    return
            elif args.run=='3':
                if args.finalstate=='mm':
                    mc = ['/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root']
                    dt = ['/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC/2022/SingleMuon_Run2022C-PromptReco-v1/mm/SingleMuon_Run2022C-PromptReco-v1_*.root']
                    dt+= ['/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC/2022/Muon_Run2022C-PromptReco-v1/mm/Muon_Run2022C-PromptReco-v1_*.root']
                    #dt+= ['/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC/2022/Muon_Run2022D-PromptReco-v1/mm/Muon_Run2022D-PromptReco-v1_*.root']
                    #dt+= ['/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC/2022/Muon_Run2022D-PromptReco-v2/mm/Muon_Run2022D-PromptReco-v2_*.root']
                elif args.finalstate=='ee':
                    mc = ['/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/ee/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root']
                    dt = ['/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf/2022/EGamma_Run2022C-PromptReco-v1/ee/EGamma_Run2022C-PromptReco-v1_*.root']
                    dt+= ['/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf/2022/EGamma_Run2022D-PromptReco-v1/ee/EGamma_Run2022D-PromptReco-v1_*.root']
                    dt+= ['/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf/2022/EGamma_Run2022D-PromptReco-v2/ee/EGamma_Run2022D-PromptReco-v2_*.root']
                else:
                    print('Only Files with two leptons of the same kind in the finalstate should be considered for corrections. Aborting...')
                    return         
            else:
                print('Run specified is neither 2 nor 3. Aborting...')
                return
            return mc, dt
        
        # paths with corrected ntuples
        if args.corr:
            if args.run=='2':
                if args.finalstate=='mm':
                    mc = ['/ceph/moh/CROWN_samples/EarlyRun3_V12/friends/mu_corr_{v}/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/mm/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_*.root'.format(v=args.version)]
                    dt = ['/ceph/moh/CROWN_samples/EarlyRun3_V12/friends/mu_corr_{v}/2018/SingleMuon_Run2018{run}-UL2018/mm/SingleMuon_Run2018{run}-UL2018_*.root'.format(run=r, v=args.version) for r in ['A', 'B', 'C', 'D']]
                elif args.finalstate=='ee':
                    mc = ['/ceph/moh/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/ee/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_*.root'.format(v=args.version)]
                    dt = [
                        '/ceph/moh/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018{run}-UL2018/ee/EGamma_Run2018{run}-UL2018_*.root'.format(run=r, v=args.version) for r in ['A', 'B', 'C', 'D']
                        ]
                    """
                    dt.append(['/ceph/moh/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i=i) for i in range(113)])
                    dt.append(['/ceph/moh/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i='113_ee')])
                    dt.append(['/ceph/moh/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i=i) for i in range(114, 117)])
                    dt.append(['/ceph/moh/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i='117_ee')])
                    dt.append(['/ceph/moh/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i=i) for i in range(118, 173)])
                    """
                else:
                    print('Only Files with two leptons of the same kind in the finalstate should be considered for corrections. Aborting...')
                    return      
            elif args.run=='3':
                res_factor_str = str(args.res_factor)
                res_factor_str = res_factor_str.replace('.', 'p')
                if args.finalstate=='mm':
                    mc = ['/ceph/jdriesch/CROWN_samples/Run3V03/ntuples_xsec_sf_lep_corr_{v}_x{x}/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root'.format(v=args.version, x='opt')]
                    dt = ['/ceph/jdriesch/CROWN_samples/Run3V03/ntuples_xsec_sf_lep_corr_{v}_x{x}/2022/SingleMuon_Run2022C-PromptReco-v1/mm/SingleMuon_Run2022C-PromptReco-v1_*.root'.format(v=args.version, x='opt')]
                    dt+= ['/ceph/jdriesch/CROWN_samples/Run3V03/ntuples_xsec_sf_lep_corr_{v}_x{x}/2022/Muon_Run2022C-PromptReco-v1/mm/Muon_Run2022C-PromptReco-v1_*.root'.format(v=args.version, x='opt')]
                    #dt+= ['/ceph/moh/CROWN_samples/Run3V02/friends/{muel}{corr}_{v}_x{x}/2022/Muon_Run2022D-PromptReco-v1/mm/Muon_Run2022D-PromptReco-v1_*.root'.format(muel='lep', corr='_corr', v=args.version, x=res_factor_str)]
                    #dt+= ['/ceph/moh/CROWN_samples/Run3V02/friends/{muel}{corr}_{v}_x{x}/2022/Muon_Run2022D-PromptReco-v2/mm/Muon_Run2022D-PromptReco-v2_*.root'.format(muel='lep', corr='_corr', v=args.version, x=res_factor_str)]
                elif args.finalstate=='ee':
                    mc = ['/ceph/moh/CROWN_samples/Run3V02/friends/{muel}{corr}_{v}_x{x}/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/ee/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root'.format(muel='lep', corr='_corr', v=args.version, x=res_factor_str)]
                    dt = ['/ceph/moh/CROWN_samples/Run3V02/friends/{muel}{corr}_{v}_x{x}/2022/EGamma_Run2022C-PromptReco-v1/ee/EGamma_Run2022C-PromptReco-v1_*.root'.format(muel='lep', corr='_corr', v=args.version, x=res_factor_str)]
                    dt+= ['/ceph/moh/CROWN_samples/Run3V02/friends/{muel}{corr}_{v}_x{x}/2022/EGamma_Run2022D-PromptReco-v1/ee/EGamma_Run2022D-PromptReco-v1_*.root'.format(muel='lep', corr='_corr', v=args.version, x=res_factor_str)]
                    dt+= ['/ceph/moh/CROWN_samples/Run3V02/friends/{muel}{corr}_{v}_x{x}/2022/EGamma_Run2022D-PromptReco-v2/ee/EGamma_Run2022D-PromptReco-v2_*.root'.format(muel='lep', corr='_corr', v=args.version, x=res_factor_str)]
                else:
                    print('Only Files with two leptons of the same kind in the finalstate should be considered for corrections. Aborting...')
                    return                         
            return mc, dt #/ceph/jdriesch/CROWN_samples/Run3V03/ntuples_xsec_sf_lep_corr_02_xopt/2022/SingleMuon_Run2022C-PromptReco-v1/mm/

    # get correction file's paths
    elif mode == 1:
        if 'mm' in args.finalstate:
            fs = 'mm'
        else:
            fs = 'ee'
        path = 'correction_files/Run{run}/{fs}/{res}mz_{mcdt}_{bin1}_{bin2}_{v}{corr}.txt'.format(run=args.run, fs=fs, bin1=kw[0], bin2=kw[1], corr=kw[2], res = kw[3], mcdt = kw[4], v=args.version)
        print("take correction file: ", path)
        return path


def calc_m(rdf, corr = ""):
    rdf = rdf.Define("lv_1", "ROOT::Math::PtEtaPhiMVector p(pt_1"+corr+", eta_1, phi_1, mass_1); return p")
    rdf = rdf.Define("lv_2", "ROOT::Math::PtEtaPhiMVector p(pt_2"+corr+", eta_2, phi_2, mass_2); return p")
    
    rdf = rdf.Define("dimuon", "lv_1 + lv_2")
    rdf = rdf.Define("pt_vis"+corr, "dimuon.Pt()")
    rdf = rdf.Define("m_vis"+corr, "dimuon.M()")
    return rdf


def plot_binned(rdf, mcdt, filters, corr='', weight='1.'):
    rdf_filtered = rdf.Filter(filters)
    #print(rdf_filtered.Count().GetValue())

    rdf_filtered_weight = rdf_filtered.Define("weight", weight)

    hist = rdf_filtered_weight.Histo1D((mcdt, 'm_vis'+corr, 200, 50, 130), 'm_vis'+corr, "weight")
    num_entries = hist.GetEntries()
    #hist.Scale(1./hist.Integral())
    hist.SetXTitle("m_vis{} (GeV)".format(corr))
    hist.SetYTitle("a.u.")

    return hist, num_entries


def make_hists(bins, args):
    bin1, bin2 = list(bins.keys())[0], list(bins.keys())[1]
    nbins1, nbins2 = len(bins[bin1])-1, len(bins[bin2])-1
    n_events = np.zeros((nbins1, nbins2))
    
    corr, corr1, corr2 = '', '', ''

    mc, dt = get_paths(args, mode=0)
    print(mc)

    ch_mc = ROOT.TChain("ntuple")
    ch_dt = ROOT.TChain("ntuple")

    for p in mc:
        ch_mc.Add(p)

    for p in dt:
        ch_dt.Add(p)

    rdf_mc = ROOT.RDataFrame(ch_mc)
    rdf_dt = ROOT.RDataFrame(ch_dt)

    event_selection = "(q_1*q_2 < 0)"
    if args.finalstate == "mm":
        event_selection += " && (trg_single_mu24_1 || trg_single_mu24_2)"
    elif args.finalstate == "ee":
        event_selection += " && (trg_single_ele27_1 || trg_single_ele27_2)"

    if int(args.run_i) > 0:
        event_selection += f" && (run >= {int(args.run_i)} || run == 1)"
    if int(args.run_f) > 0:
        event_selection += f" && (run <= {int(args.run_f)} || run == 1)"

    rdf_mc = rdf_mc.Filter(event_selection)
    rdf_dt = rdf_dt.Filter(event_selection)    

    if args.corr:
        if bin1=='pt':
            corr1 = '_corr'
            corr = '_corr'
        if bin2=='pt':
            corr2='_corr'
            corr = '_corr'


    outdir = 'hists/Run{run}/{fs}/{bin1}_{bin2}{corr}_{v}/'.format(run=args.run, fs=args.finalstate, bin1=bin1, bin2=bin2, corr=corr, v=args.version)
    if not utils.usedir(outdir, args.overwrite):
        return


    with open(outdir+"bins.json", 'w') as fp:
        json.dump(bins, fp)

    for i in range(nbins1):
        for j in tqdm(range(nbins2)):
            file0 = ROOT.TFile(outdir+'{bin1}_{}_{bin2}_{}.root'.format(i,j, bin1=bin1, bin2=bin2), 'RECREATE')

            bin1_l, bin1_r = bins[bin1][i], bins[bin1][i+1]
            bin2_l, bin2_r = bins[bin2][j], bins[bin2][j+1]

            # filter for a muon in event in corresponding bins 
            filter_template = "({bin}_{n}{corr} > {bin_l} && {bin}_{n}{corr} < {binr})"

            filter1a = filter_template.format(bin=bin1, n=1, bin_l=bin1_l, binr = bin1_r, corr=corr1)
            filter1b = filter_template.format(bin=bin2, n=1, bin_l=bin2_l, binr = bin2_r, corr=corr2)
            
            filter2a = filter_template.format(bin=bin1, n=2, bin_l=bin1_l, binr = bin1_r, corr=corr1)
            filter2b = filter_template.format(bin=bin2, n=2, bin_l=bin2_l, binr = bin2_r, corr=corr2)

            filters = "({} && {}) || ({} && {})".format(filter1a, filter1b, filter2a, filter2b)

            # plot di muon mass in filtered events
            hist_mc, num_mc = plot_binned(rdf_mc, "mc", filters, corr, "genweight*sumwWeight*crossSectionPerEventWeight*sf_trk*sf_sta*sf_id*sf_iso*sf_trg")
            hist_dt, num_dt = plot_binned(rdf_dt, "data", filters, corr)
            n_events[i][j] = num_dt
            if args.info:
                print("###############################################################################")
                print("{} < {} < {} ,  {} < {} < {}".format(bins[bin1][i], bin1, bins[bin1][i+1], bins[bin2][j], bin2, bins[bin2][j+1]))
                print("mc entries: {} ,  dt entries: {}".format(num_mc, num_dt))
                print("")
            hist_mc.Write()
            hist_dt.Write()
            file0.Close()
    np.savetxt(outdir+'dt_events.txt', n_events)


def get_corrections(args, bins):
    if args.corr:
        corr='_corr'
    else:
        corr=''
    bin1, bin2 = list(bins.keys())[0], list(bins.keys())[1]
    nbins1, nbins2 = len(bins[bin1])-1, len(bins[bin2])-1
    mz_mc, mz_dt, res_mz_mc, res_mz_dt = np.zeros((nbins1, nbins2)), np.zeros((nbins1, nbins2)), np.zeros((nbins1, nbins2)), np.zeros((nbins1, nbins2))
    k = 0
    
    outdir = 'plots/Run{run}/{fs}/{bin1}_{bin2}{corr}_{v}/'.format(run=args.run, fs=args.finalstate, bin1=bin1, bin2=bin2, corr=corr, v=args.version)
    if not utils.usedir(outdir, args.overwrite):
        return
    print("Calculating scale and resolution corrections")
    for i in range(len(bins[bin1])-1):
        for j in range(len(bins[bin2])-1):
            
            fit_results = utils.roofit_mass('hists/Run{run}/{fs}/{bin1}_{bin2}{corr}_{v}/{bin1}_{}_{bin2}_{}.root'.format(i, j, run=args.run, fs=args.finalstate, bin1=bin1, bin2=bin2, corr=corr, v=args.version), fitf='bwxcb', plot=outdir+'{}_{}'.format(i, j))
            
            # make sure that chi2 is not through the roof
            c2 = 0
            while fit_results['mc'][2]>70 or fit_results['dt'][2]>100:
                fit_results = utils.roofit_mass('hists/Run{run}/{fs}/{bin1}_{bin2}{corr}_{v}/{bin1}_{}_{bin2}_{}.root'.format(i, j, run=args.run, fs=args.finalstate, bin1=bin1, bin2=bin2, corr=corr, v=args.version), fitf='bwxcb', plot=outdir+'{}_{}'.format(i, j))
                c2+=1
                if c2 > 10:
                    print('Fit did not work for 10 consecutive times: Chi2/dof>70')
                    return
            mz_mc[i][j] = fit_results['mc'][0]
            mz_dt[i][j] = fit_results['dt'][0]
            res_mz_mc[i][j] = fit_results['mc'][1] # account for sigma_l, sigma_r
            res_mz_dt[i][j] = fit_results['dt'][1]
            k+=1

        print('done {} of {}'.format(k, nbins1*nbins2))  
    
    if not utils.usedir('correction_files/Run{}/{}/'.format(args.run, args.finalstate), args.overwrite):
        return
        
    np.savetxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, '', 'mc']), mz_mc)
    np.savetxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, '', 'dt']), mz_dt)

    np.savetxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, 'res', 'mc']), res_mz_mc)
    np.savetxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, 'res', 'dt']), res_mz_dt)


    # make 2d plots
    ratio_mu = (91.1876+mz_dt)/(91.1876+mz_mc)
    ratio_sigma = res_mz_dt/res_mz_mc * ratio_mu
    labels1, labels2 = np.zeros(nbins1), np.zeros(nbins2)
    
    for i in range(nbins1):
        labels1[i] = .5*(bins[bin1][i] + bins[bin1][i+1])

    for j in range(nbins2):
        labels2[j] = .5*(bins[bin2][j] + bins[bin2][j+1])

    if args.finalstate == 'ee':
        mean_min, mean_max = 0.98, 1.02
        std_min, std_max = 0.8, 1.2
    elif args.finalstate == 'mm':
        mean_min, mean_max = 0.998, 1.002
        std_min, std_max = 0.9, 1.1
    else:
        print("This finalstate should not be used to obtain lepton corrections")
        return
    utils.plot2d(ratio_mu, '{}mu{}.pdf'.format(outdir,corr), r'$\frac{M_Z (\mathrm{data})}{M_Z (\mathrm{mc})}$', r'$p_\mathrm{T}$', bins[bin2], r'$\eta$', bins[bin1], mean_min, mean_max, np.around(labels2, 1), labels1)
    utils.plot2d(ratio_sigma, '{}sigma{}.pdf'.format(outdir,corr), r'$\frac{\mathrm{res}(M_Z(\mathrm{data}))}{\mathrm{res}(M_Z (\mathrm{mc}))}$', r'$p_\mathrm{T}$', bins[bin2], r'$\eta$', bins[bin1], std_min, std_max, np.around(labels2,1), labels1)


ROOT.gInterpreter.Declare("""
        float gaus(){
            return gRandom->Gaus(0,1);
        }
        """)


def apply_corrections(args, bins):
    bin1, bin2 = list(bins.keys())[0], list(bins.keys())[1]
    nbins1, nbins2 = len(bins[bin1])-1, len(bins[bin2])-1

    if args.corr:
        corr='_corr'
    else:
        corr=''

    x = np.loadtxt('correction_files/Run{}/{}/sf_extra_mc_{}.txt'.format(args.run, args.finalstate, args.version))
    # Load correction files
    mz_mc, mz_dt = np.loadtxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, '', 'mc'])), np.loadtxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, '', 'dt']))
    pt_sf = (91.1876+mz_mc) / (91.1876+mz_dt)
    mz_res_mc, mz_res_dt = np.loadtxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, 'res', 'mc'])), np.loadtxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, 'res', 'dt'])) * pt_sf

    period = args.inpath.split('/')[4]
    process = args.inpath.split('/')[-2]
    year = args.inpath.split('/')[-3]

    inpath = args.inpath + "/" + args.finalstate + "/" + process
    
    if args.test:
        outdir = '/ceph/jdriesch/CROWN_samples/test/'
    else:
        helpdir = '/ceph/jdriesch/CROWN_samples/{period}/ntuples_xsec_sf_{muel}{corr}_{v}_x{x}/{year}/{proc}/{f}/'
        if args.finalize:
            #res_factor_str = f"{float(args.res_factor):.2f}"
            #res_factor_str = res_factor_str.replace('.', 'p')
            res_factor_str = 'opt'
            outdir = helpdir.format(period=period, muel='lep', corr='_corr', v=args.version, x=res_factor_str, year=year, proc=process, f=args.finalstate)
        else:
            if 'mm' in args.finalstate:
                muel = 'mu_corr'
            else:
                muel = 'el_corr'
            outdir = helpdir.format(period=period, muel=muel, corr=corr, v=args.version, year=year, proc=process, f=args.finalstate)


    print("Input samples taken from: {}".format(inpath))

    # for other processes than data: apply resolution correction
    if (not ("Muon" in process)) and (not ("EGamma" in process)):
        print("MC process: resolution corrections will be applied")
        data=False
    else:
        print("Data process: scale corrections will be applied")
        data=True

    if args.finalstate == "mm" or args.finalstate == 'ee':
        print("Di-lepton finalstate: both leptons will be corrected")
        dilepton=True
    else:
        print("Single-lepton finalstate: one lepton will be corrected")
        dilepton=False

    print("Output samples will be saved in: {}".format(outdir))


    # check if outdir exists and create directory if not 
    if not utils.usedir(outdir, args.overwrite):
        return

    # get list of all files in the input directory
    files = glob.glob(inpath + "_*.root")


    for f in tqdm(files):
        outfile = outdir + f.split("/")[-1]

        rdf = ROOT.RDataFrame('ntuple', f)
        original_cols = [str(col) for col in rdf.GetColumnNames()]

        rdf = rdf.Define("pt_1_corr", "pt_1")
        if dilepton:
            rdf = rdf.Define("pt_2_corr", "pt_2")

        for j in range(nbins1): 
            for k in range(nbins2): 
                bin1_l, bin1_r = bins[bin1][j], bins[bin1][j+1]
                bin2_l, bin2_r = bins[bin2][k], bins[bin2][k+1]

                # filter for a muon in event in corresponding bins
                filter_template = "({bin}_{n} > {bin_l} && {bin}_{n} < {binr})"

                filter1a = filter_template.format(bin=bin1, n=1, bin_l=bin1_l, binr = bin1_r)
                filter1b = filter_template.format(bin=bin2, n=1, bin_l=bin2_l, binr = bin2_r)

                if mz_res_mc[j][k] > mz_res_dt[j][k]:
                    res_sf = 1
                else:
                    res_sf = mz_res_dt[j][k] / mz_res_mc[j][k]


                if data:
                    rdf = rdf.Redefine("pt_1_corr", "double p; if ({} && {}) p=pt_1 * (91.1876 + {})/(91.1876 + {}); else p=pt_1_corr; return p;".format(filter1a, filter1b, str(mz_mc[j][k]), str(mz_dt[j][k])))
                else:
                    rdf = rdf.Redefine("pt_1_corr", "double p; if ({} && {}) p=pt_1 * (1 + {}/91.1876*sqrt({}-1)*(float)(gaus())); else p=pt_1_corr; return p;".format(filter1a, filter1b, str(x[j][k]*mz_res_mc[j][k]), str((res_sf)**2)))


                if dilepton:
                    filter2a = filter_template.format(bin=bin1, n=2, bin_l=bin1_l, binr = bin1_r)
                    filter2b = filter_template.format(bin=bin2, n=2, bin_l=bin2_l, binr = bin2_r)
                    
                    if data:
                        rdf = rdf.Redefine("pt_2_corr", "double p; if ({} && {}) p=pt_2 * (91.1876 + {})/(91.1876 + {}); else p=pt_2_corr; return p;".format(filter2a, filter2b,  str(mz_mc[j][k]), str(mz_dt[j][k])))
                    else:
                        rdf = rdf.Redefine("pt_2_corr", "double p; if ({} && {}) p=pt_2 * (1 + {}/91.1876*sqrt({}-1)*(float)(gaus())); else p=pt_2_corr; return p;".format(filter2a, filter2b, str(x[j][k]*mz_res_mc[j][k]), str((res_sf)**2)))

        if dilepton:
            rdf = calc_m(rdf, "_corr")
            if args.finalize:
                quants =  ["pt_1_corr", "pt_2_corr", "m_vis_corr", "pt_vis_corr"]
            else:
                quants =  ["pt_1_corr", "pt_2_corr", "eta_1", "eta_2", "phi_1", "phi_2", "mass_1", "mass_2", "pt_1", "pt_2", "m_vis", "m_vis_corr", "q_1", "q_2", "pt_vis", "pt_vis_corr"]
                if args.finalstate=='mm':
                    quants+=["trg_single_mu24_1", "trg_single_mu24_2"]
                else:
                    quants+=["trg_single_ele27_1", "trg_single_ele27_2"]

        else:
            if args.finalize:
                quants =  ["pt_1_corr"]
            else:
                quants = ["pt_1_corr", "eta_1", "phi_1", "mass_1", "pt_1", "q_1"]

        # add recoil variables
        if not ("WtoLNu" in outfile):
            bosonphi = "phi_vis_c"
            bosonpt = "pt_vis_c"
        else:
            bosonphi = "genbosonphi"
            bosonpt = "genbosonpt"

        if args.finalstate == "ee" or args.finalstate == "mm":
            rdf = rdf.Define("pt_vis_c_x", "pt_1_corr*cos(phi_1) + pt_2_corr*cos(phi_2)")
            rdf = rdf.Define("pt_vis_c_y", "pt_1_corr*sin(phi_1) + pt_2_corr*sin(phi_2)")
            rdf = rdf.Define("pt_vis_c", "sqrt(pt_vis_c_x*pt_vis_c_x + pt_vis_c_y*pt_vis_c_y)")
            rdf = rdf.Define("phi_vis_c", "atan2(pt_vis_c_y, pt_vis_c_x)")
        else:
            rdf = rdf.Define("pt_vis_c", "pt_1_corr")
            rdf = rdf.Define("phi_vis_c", "phi_1")

        rdf = rdf.Define("bosonpt", f"{bosonpt}")
        rdf = rdf.Define("bosonphi", f"{bosonphi}")

        rdf = rdf.Define("uPx", "met_uncorrected*cos(metphi_uncorrected) + pt_vis_c*cos(phi_vis_c)")
        rdf = rdf.Define("uPy", "met_uncorrected*sin(metphi_uncorrected) + pt_vis_c*sin(phi_vis_c)")

        rdf = rdf.Define("uP1_uncorrected", "- (uPx*cos("+bosonphi+") + uPy*sin("+bosonphi+"))")
        rdf = rdf.Define("uP2_uncorrected", "uPx*sin("+bosonphi+") - uPy*cos("+bosonphi+")")

        rdf = rdf.Define("pfuPx", "pfmet_uncorrected*cos(pfmetphi_uncorrected) + pt_vis_c*cos(phi_vis_c)")
        rdf = rdf.Define("pfuPy", "pfmet_uncorrected*sin(pfmetphi_uncorrected) + pt_vis_c*sin(phi_vis_c)")

        rdf = rdf.Define("pfuP1_uncorrected", "- (pfuPx*cos("+bosonphi+") + pfuPy*sin("+bosonphi+"))")
        rdf = rdf.Define("pfuP2_uncorrected", "pfuPx*sin("+bosonphi+") - pfuPy*cos("+bosonphi+")")

        met_cols = ["uP1_uncorrected", "uP2_uncorrected", "pfuP1_uncorrected", "pfuP2_uncorrected", "pt_vis_c", "phi_vis_c", "bosonpt", "bosonphi"]

        #print(original_cols + quants + met_cols)
        rdf.Snapshot("ntuple", outfile, original_cols + quants + met_cols)

    print("Great success!")



def plot(args, bins):
    corr=''
    if args.corr:
        corr='_corr'
    bin1, bin2 = list(bins.keys())[0], list(bins.keys())[1]
    nbins1, nbins2 = len(bins[bin1])-1, len(bins[bin2])-1
    outdir = 'plots/Run{run}/{fs}/ratio_plots{corr}_{bin1}_{bin2}_{v}/'.format(run=args.run, fs=args.finalstate, corr=corr, bin1=bin1, bin2=bin2, v=args.version)
    if not utils.usedir(outdir, args.overwrite):
        return

    chi2 = np.zeros((nbins1, nbins2))
    for i in range(nbins1):
        for j in range(nbins2):
            filename = 'hists/Run{run}/{fs}/{bin1}_{bin2}{corr}_{v}/{bin1}_{}_{bin2}_{}.root'.format(i, j, run=args.run, fs=args.finalstate, bin1=bin1, bin2=bin2, corr=corr, v=args.version)
            file0 = ROOT.TFile(filename)
            hist_mc, hist_dt = file0.Get('mc'), file0.Get('data')
            hist_mc.Scale(1./hist_mc.Integral())
            n_dt = round(hist_dt.Integral())
            hist_dt.Scale(1./hist_dt.Integral())
            hists = {'mc': hist_mc, 'dt': hist_dt}
            outfile = outdir+'{b1}_{b2}.pdf'.format(b1=i, b2=j)
            text = ['lepton 1 or 2:', '{} < #eta < {}'.format(bins[bin1][i], bins[bin1][i+1]), '{} GeV < pT < {} GeV'.format(bins[bin2][j], bins[bin2][j+1])]
            chi2[i][j] = utils.plot_ratio(plots=hists, rcolors={'mc': ROOT.kBlue, 'dt': ROOT.kBlack}, title='ratio plot', outfile=outfile, evts=n_dt, text=text)
    np.savetxt(outdir+'chi2.txt', chi2)
    

if __name__=='__main__':
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    args = parse_args()

    if 'mm' in args.finalstate:
        bins = {'eta': [-2.4, -1.2, 0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]}
    else:
        bins = {'eta': [-2.5, -1.49, 0., 1.49, 2.5], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]}
    
    if args.histogram:
        make_hists(bins, args)

    if args.calculation:
        get_corrections(args, bins)
        #get_pt_rel_res(bins)

    if args.evaluation:
        apply_corrections(args, bins)

    if args.plot:
        plot(args, bins)
