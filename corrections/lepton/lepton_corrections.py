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
                    mc = ['/ceph/moh/CROWN_samples/Run3V01/ntuples/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root']
                    dt = ['/ceph/moh/CROWN_samples/Run3V01/ntuples/2022/SingleMuon_Run2022C-PromptReco-v1/mm/SingleMuon_Run2022C-PromptReco-v1_*.root'] 
                elif args.finalstate=='ee':
                    mc = ['/ceph/moh/CROWN_samples/Run3V01/ntuples/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/ee/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root']
                    dt = ['/ceph/moh/CROWN_samples/Run3V01/ntuples/2022/EGamma_Run2022C-PromptReco-v1/ee/EGamma_Run2022C-PromptReco-v1_*.root']
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
                    mc = ['/ceph/jdriesch/CROWN_samples/EarlyRun3_V12/friends/mu_corr_{v}/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/mm/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_*.root'.format(v=args.version)]
                    dt = ['/ceph/jdriesch/CROWN_samples/EarlyRun3_V12/friends/mu_corr_{v}/2018/SingleMuon_Run2018{run}-UL2018/mm/SingleMuon_Run2018{run}-UL2018_*.root'.format(run=r, v=args.version) for r in ['A', 'B', 'C', 'D']]
                elif args.finalstate=='ee':
                    mc = ['/ceph/jdriesch/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/ee/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X_*.root'.format(v=args.version)]
                    dt = [
                        '/ceph/jdriesch/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018{run}-UL2018/ee/EGamma_Run2018{run}-UL2018_*.root'.format(run=r, v=args.version) for r in ['A', 'B', 'C', 'D']
                        ]
                    """
                    dt.append(['/ceph/jdriesch/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i=i) for i in range(113)])
                    dt.append(['/ceph/jdriesch/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i='113_ee')])
                    dt.append(['/ceph/jdriesch/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i=i) for i in range(114, 117)])
                    dt.append(['/ceph/jdriesch/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i='117_ee')])
                    dt.append(['/ceph/jdriesch/CROWN_samples/EarlyRun3_V12/friends/el_corr_{v}/2018/EGamma_Run2018D-UL2018/ee/EGamma_Run2018D-UL2018_{i}.root'.format(v=args.version, i=i) for i in range(118, 173)])
                    """
                else:
                    print('Only Files with two leptons of the same kind in the finalstate should be considered for corrections. Aborting...')
                    return      
            elif args.run=='3':
                if args.finalstate=='mm':
                    mc = ['/ceph/jdriesch/CROWN_samples/Run3V01/friends/mu_corr_{v}/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root'.format(v=args.version)]
                    dt = ['/ceph/jdriesch/CROWN_samples/Run3V01/friends/mu_corr_{v}/2022/SingleMuon_Run2022C-PromptReco-v1/mm/SingleMuon_Run2022C-PromptReco-v1_*.root'.format(v=args.version)]
                elif args.finalstate=='ee':
                    mc = ['/ceph/jdriesch/CROWN_samples/Run3V01/friends/el_corr_{v}/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/ee/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root'.format(v=args.version)]
                    dt = ['/ceph/jdriesch/CROWN_samples/Run3V01/friends/el_corr_{v}/2022/EGamma_Run2022C-PromptReco-v1/ee/EGamma_Run2022C-PromptReco-v1_*.root'.format(v=args.version)]
                else:
                    print('Only Files with two leptons of the same kind in the finalstate should be considered for corrections. Aborting...')
                    return                         
            return mc, dt

    # get correction file's paths
    elif mode == 1:
        if 'mm' in args.finalstate:
            fs = 'mm'
        else:
            fs = 'ee'
        path = 'correction_files/Run{run}/{fs}/{res}mz_{mcdt}_{bin1}_{bin2}_{v}{corr}.txt'.format(run=args.run, fs=fs, bin1=kw[0], bin2=kw[1], corr=kw[2], res = kw[3], mcdt = kw[4], v=args.version)
        return path


def calc_m(rdf, corr = ""):
    rdf = rdf.Define("lv_1", "ROOT::Math::PtEtaPhiMVector p(pt_1"+corr+", eta_1, phi_1, mass_1); return p")
    rdf = rdf.Define("lv_2", "ROOT::Math::PtEtaPhiMVector p(pt_2"+corr+", eta_2, phi_2, mass_2); return p")
    
    rdf = rdf.Define("dimuon", "lv_1 + lv_2")
    rdf = rdf.Define("m_vis"+corr, "dimuon.M()")
    return rdf


def plot_binned(rdf, mcdt, filters, corr=''):
    rdf_filtered = rdf.Filter(filters)
    #print(rdf_filtered.Count().GetValue())

    hist = rdf_filtered.Histo1D((mcdt, 'm_vis'+corr, 200, 50, 130), 'm_vis'+corr)
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

    rdf_mc = ROOT.RDataFrame("ntuple", mc)
    rdf_dt = ROOT.RDataFrame("ntuple", dt)

    if args.corr:
        rdf_mc = calc_m(rdf_mc, '_corr')
        rdf_dt = calc_m(rdf_dt, '_corr')
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
            hist_mc, num_mc = plot_binned(rdf_mc, "mc", filters, corr)
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
    ratio_mu = mz_dt/mz_mc
    ratio_sigma = res_mz_dt/res_mz_mc * ratio_mu
    labels1, labels2 = np.zeros(nbins1), np.zeros(nbins2)
    
    for i in range(nbins1):
        labels1[i] = .5*(bins[bin1][i] + bins[bin1][i+1])

    for j in range(nbins2):
        labels2[j] = .5*(bins[bin2][j] + bins[bin2][j+1])

    
    utils.plot2d(ratio_mu, '{}mu{}.pdf'.format(outdir,corr), r'$\frac{M_Z (\mathrm{data})}{M_Z (\mathrm{mc})}$', r'$p_\mathrm{T}$', bins[bin2], r'$\eta$', bins[bin1], 0.998, 1.002, np.around(labels2, 1), labels1)
    utils.plot2d(ratio_sigma, '{}sigma{}.pdf'.format(outdir,corr), r'$\frac{\mathrm{res}(M_Z(\mathrm{data}))}{\mathrm{res}(M_Z (\mathrm{mc}))}$', r'$p_\mathrm{T}$', bins[bin2], r'$\eta$', bins[bin1], 0.9, 1.1, np.around(labels2,1), labels1)


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

    # Load correction files
    #fname = get_paths(args, mode=1, kw=[bin1, bin2, corr])
    mz_mc, mz_dt = np.loadtxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, '', 'mc'])), np.loadtxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, '', 'dt']))
    pt_sf = mz_mc / mz_dt
    mz_res_mc, mz_res_dt = np.loadtxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, 'res', 'mc'])), np.loadtxt(get_paths(args, mode=1, kw=[bin1, bin2, corr, 'res', 'dt'])) * pt_sf

    period = args.inpath.split('/')[4]
    process = args.inpath.split('/')[-2]
    year = args.inpath.split('/')[-3]

    inpath = args.inpath + "/" + args.finalstate + "/" + process
    
    if args.test:
        outdir = '/ceph/jdriesch/CROWN_samples/test/'
    else:
        helpdir = '/ceph/jdriesch/CROWN_samples/{period}/friends/{muel}{corr}_{v}/{year}/{proc}/{f}/'
        if args.finalize:
            outdir = helpdir.format(period=period, muel='lep', corr='', v='corr', year=year, proc=process, f=args.finalstate)
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
                    rdf = rdf.Redefine("pt_1_corr", "double p; if ({} && {}) p=pt_1 * {}; else p=pt_1_corr; return p;".format(filter1a, filter1b, str(pt_sf[j][k])))
                else:
                    rdf = rdf.Redefine("pt_1_corr", "double p; if ({} && {}) p=pt_1 * (1 + {}/91*sqrt({}-1)*(float)(gaus())); else p=pt_1_corr; return p;".format(filter1a, filter1b, str(mz_res_mc[j][k]), str((res_sf)**2))) # TODO: evaluate if division by 2 needed


                if dilepton:
                    filter2a = filter_template.format(bin=bin1, n=2, bin_l=bin1_l, binr = bin1_r)
                    filter2b = filter_template.format(bin=bin2, n=2, bin_l=bin2_l, binr = bin2_r)
                    
                    if data:
                        rdf = rdf.Redefine("pt_2_corr", "double p; if ({} && {}) p=pt_2 * {}; else p=pt_2_corr; return p;".format(filter2a, filter2b, str(pt_sf[j][k])))
                    else:
                        rdf = rdf.Redefine("pt_2_corr", "double p; if ({} && {}) p=pt_2* (1 + {}/91*sqrt({}-1)*(float)(gaus())); else p=pt_2_corr; return p;".format(filter2a, filter2b, str(mz_res_mc[j][k]), str((res_sf)**2)))
        
        if dilepton:
            rdf = calc_m(rdf, "_corr")
            if args.finalize:
                quants =  ["pt_1_corr", "pt_2_corr", "m_vis_corr"]
            else:
                quants =  ["pt_1_corr", "pt_2_corr", "eta_1", "eta_2", "phi_1", "phi_2", "mass_1", "mass_2", "pt_1", "pt_2", "m_vis_corr"]
        
        else:
            if args.finalize:
                quants =  ["pt_1_corr"]
            else:
                quants = ["pt_1_corr", "eta_1", "phi_1", "mass_1", "pt_1"]

        rdf.Snapshot("ntuple", outfile, quants)

    print("Great success!")



def plot(args, bins):
    corr=''
    if args.corr:
        corr='_corr'
    bin1, bin2 = list(bins.keys())[0], list(bins.keys())[1]
    nbins1, nbins2 = len(bins[bin1])-1, len(bins[bin2])-1
    if not utils.usedir('plots/Run{run}/{fs}/ratio_plots{corr}_{bin1}_{bin2}_{v}/'.format(run=args.run, fs=args.finalstate, corr=corr, bin1=bin1, bin2=bin2, v=args.version), args.overwrite):
        return
    # Load files
    for i in range(nbins1):
        for j in range(nbins2):
            filename = 'hists/Run{run}/{fs}/{bin1}_{bin2}{corr}_{v}/{bin1}_{}_{bin2}_{}.root'.format(i, j, run=args.run, fs=args.finalstate, bin1=bin1, bin2=bin2, corr=corr, v=args.version)
            file0 = ROOT.TFile(filename)
            hist_mc, hist_dt = file0.Get('mc'), file0.Get('data')
            hist_mc.Scale(1./hist_mc.Integral())
            hist_dt.Scale(1./hist_dt.Integral())
            hists = {'mc': hist_mc, 'dt': hist_dt}
            outfile = 'plots/Run{run}/{fs}/ratio_plots{corr}_{bin1}_{bin2}_{v}/{b1}_{b2}.pdf'.format(run=args.run, bin1=bin1, bin2=bin2, b1=i, b2=j, fs=args.finalstate, corr=corr, v=args.version)
            utils.plot_ratio(plots=hists, rcolors={'mc': ROOT.kBlue, 'dt': ROOT.kBlack}, title='ratio plot', outfile=outfile) #, ratio=True)


if __name__=='__main__':
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    args = parse_args()

    bins = {'eta': [-2.4, -1.2, 0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 1500]}
    
    if args.histogram:
        make_hists(bins, args)

    if args.calculation:
        get_corrections(args, bins)
        #get_pt_rel_res(bins)

    if args.evaluation:
        apply_corrections(args, bins)

    if args.plot:
        plot(args, bins)
