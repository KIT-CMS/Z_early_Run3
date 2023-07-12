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
from multiprocessing import Pool, RLock


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
            Lepton correction program. To obtain corrections, please proceed 
            as follows:\n 
            1. Produce histograms by setting option -H\n 
            2. Fit Z mass peak position and width via -F\n 
            3. Produce smeared ntuples and histograms via -S\n 
            4. Fit histograms with extra smearing factor via -C\n
            5. Plot fit results and calculate extra smearing via -P
            """)
    parser.add_argument(
        '-H',
        '--histogram',
        action='store_true',
        default=False,
        help='Create binned histograms'
    )
    parser.add_argument(
        '-F',
        '--fit',
        action='store_true',
        default=False,
        help='Mass fit to histograms'
    )
    parser.add_argument(
        '-S',
        '--smear',
        action='store_true',
        default=False,
        help='make ntuples for smearing factor'
    )
    parser.add_argument(
        '-C',
        '--calculation',
        action='store_true',
        default=False,
        help='set to calculate extra smearing'
    )
    parser.add_argument(
        '-P',
        '--plot',
        action='store_true',
        default=False,
        help='set to calculate extra smearing'
    )
    parser.add_argument(
        '--corr',
        action='store_true',
        default=False
    )
    parser.add_argument('--run_i', default=355862)
    parser.add_argument('--run_f', default=357482)
    args = parser.parse_args()
    return args



def calc_m(rdf, corr = ""):
    """
    function for calculation of pt and mass of a dilepton system

    (ROOT.RDataFrame) rdf: dataframe
    (str) corr: input flag for corrected values
    """
    
    rdf = rdf.Define(
        "lv_1"+corr,
        "ROOT::Math::PtEtaPhiMVector p(pt_1"+corr+", eta_1, phi_1, mass_1);"\
        "return p"
    )
    rdf = rdf.Define(
        "lv_2"+corr,
        "ROOT::Math::PtEtaPhiMVector p(pt_2"+corr+", eta_2, phi_2, mass_2);"\
        "return p"
    )
    
    rdf = rdf.Define("dimuon"+corr, f"lv_1{corr} + lv_2{corr}")
    rdf = rdf.Define("pt_vis"+corr, f"dimuon{corr}.Pt()")
    rdf = rdf.Define("m_vis"+corr, f"dimuon{corr}.M()")
    
    return rdf


def load_tchain(rfiles, ch, corr=''):
    """
    function to load rdf

    (dict) rfiles: dictionary with channels and corresponding files
    (str) ch: current channel
    """

    tchain = ROOT.TChain("ntuple")
    sfchain = ROOT.TChain("ntuple")
    puchain = ROOT.TChain("ntuple")
    xsecchain = ROOT.TChain("ntuple")
    lepchain = ROOT.TChain("ntuple")
    for f in rfiles[ch]:
        tchain.Add(f)
        sfchain.Add(f.replace("ntuples", "friends/sf"))
        puchain.Add(f.replace("ntuples", "friends/pu"))
        xsecchain.Add(f.replace("ntuples", "friends/xsec"))
        if corr=='_corr':
            lepchain.Add(f.replace("ntuples", "friends/lepton"))
    tchain.AddFriend(sfchain)
    tchain.AddFriend(puchain)
    tchain.AddFriend(xsecchain)
    if corr=='_corr':
        tchain.AddFriend(lepchain)
        print("Added friend with lepton momentum corrections")
    
    return tchain


def select_events(rdf, ch, args):
    """
    function to perform event selection

    (ROOT.RDataFrame) rdf: root RDataFrame with data or mc
    (str) ch: curren channel
    """

    event_selection = "(q_1*q_2 < 0)"
    if 'mm' in ch: # finalstate mm
        event_selection += " && (trg_single_mu24_1 || trg_single_mu24_2)"
    elif 'ee' in ch: # finalstate ee
        event_selection += " && (trg_single_ele27_1 || trg_single_ele27_2)"

    if int(args.run_i) > 0:
        event_selection += f" && (run >= {int(args.run_i)} || run == 1)"
    if int(args.run_f) > 0:
        event_selection += f" && (run <= {int(args.run_f)} || run == 1)"

    rdf = rdf.Filter(event_selection)

    # include weights for mc
    if 'mc' in ch:
        weight = "genweight*sumwWeight*crossSectionPerEventWeight"\
                    "*puweight*sf_trk*sf_sta*sf_id*sf_iso*sf_trg"
    else:
        weight = "1"

    rdf = rdf.Define("weight", weight)

    return rdf


def make_hists(datasets, ch, bins, args, corr=''):
    """
    function for creating histograms

    (dict) datasets: dictionary containing arrays w/ sample paths
    (str) ch: current channel
    (dict) bins: dictionary with bin edges for corrections
    (parser) args: arguments from parser
    """

    # prepare dataframe
    tchain = load_tchain(datasets, ch, corr)        
    rdf = ROOT.RDataFrame(tchain)
    rdf = select_events(rdf, ch, args)
    print(rdf.Mean("puweight").GetValue())

    # prepare storing
    outdir = f'root_files/'
    os.makedirs(outdir, exist_ok=True)
    hists = []

    for e in tqdm(range(len(bins['eta'])-1)):
        for p in range(len(bins['pt'])-1):

            # filter for a muon in event in corresponding bins 
            eta_l, eta_r = bins['eta'][e], bins['eta'][e+1]
            pt_l, pt_r = bins['pt'][p], bins['pt'][p+1]

            f_eta_1 = f"(eta_1 > {eta_l} && eta_1 < {eta_r})"
            f_pt_1 = f"(pt_1{corr} > {pt_l} && pt_1{corr} < {pt_r})"
            f_eta_2 = f_eta_1.replace('eta_1', 'eta_2')
            f_pt_2 = f_pt_1.replace('pt_1', 'pt_2')
            
            rdf_filtered = rdf.Filter(
                f"({f_eta_1} && {f_pt_1}) || ({f_eta_2} && {f_pt_2})"
                )

            # make histogram
            hist = rdf_filtered.Histo1D(
                (f"eta{e}pt{p}", 'm_vis'+corr, 200, 50, 130),
                'm_vis'+corr,
                "weight"
            )
            hist.SetXTitle("m_vis{} (GeV)".format(corr))
            hist.SetYTitle("a.u.")
            hists.append(hist)

    file0 = ROOT.TFile(outdir+f"{ch}_hists{corr}.root", 'RECREATE')
    for hist in hists:
        hist.Write()
    file0.Close()


def plot_fits(n_eta, n_pt, corr=''):
    """
    function to fit mass histogram to bw x cb and save results
    
    (dict) bins: dictionary with histogram bin edges
    TODO add corrected hists, generalize for ee
    """
    mz_mc, res_mc = np.zeros((n_eta, n_pt)), np.zeros((n_eta, n_pt))
    mz_dt, res_dt = np.zeros((n_eta, n_pt)), np.zeros((n_eta, n_pt))

    hfile_mc = ROOT.TFile(f'root_files/mcmm_hists{corr}.root')
    hfile_dt = ROOT.TFile(f'root_files/2022Cmm_hists{corr}.root')

    plotdir = f"plots/unsmeared_massfits{corr}"
    os.makedirs(plotdir, exist_ok=True)

    arguments = []
    for e in range(n_eta):
        for p in range(n_pt):
            arguments.append((
                hfile_mc.Get(f"eta{e}pt{p}"),
                hfile_dt.Get(f"eta{e}pt{p}"),
                f'{plotdir}/eta{e}pt{p}',
                'bwxcb',
                False,
                e,
                p,
            ))

    pool = Pool(8, initargs=(RLock(),), initializer=tqdm.set_lock)
    for results in tqdm(pool.imap_unordered(job_wrapper_fits, arguments)):
        fit, fit_err, e, p = results
        mz_mc[e][p], res_mc[e][p] = fit['mc'][0], fit['mc'][1]
        mz_dt[e][p], res_dt[e][p] = fit['dt'][0], fit['dt'][1]

    corr_dir = f"correction_files/"
    np.savetxt(corr_dir+f"mc{corr}_response.txt", mz_mc)
    np.savetxt(corr_dir+f"mc{corr}_resolution.txt", res_mc)
    np.savetxt(corr_dir+f"dt{corr}_response.txt", mz_dt)
    np.savetxt(corr_dir+f"dt{corr}_resolution.txt", res_dt)


def make_ratio(n_eta, n_pt, corr=''):
    hfile_mc = ROOT.TFile(f'root_files/mcmm_hists{corr}.root')
    hfile_dt = ROOT.TFile(f'root_files/2022Cmm_hists{corr}.root')

    for e in range(n_eta):
        for p in range(n_pt):
            plots = {
                'mc': hfile_mc.Get(f'eta{e}pt{p}'),
                'dt': hfile_dt.Get(f'eta{e}pt{p}')
            }
            chi2 = utils.plot_ratio(
                plots=plots, 
                title='ratio plot', 
                outfile=f'plots/unsmeared_massfits{corr}/eta{e}pt{p}_ratio', 
                xrange=[80,102]
            )


def make2d(bins, corr=''):
    """
    function for 2d plot of resolution and response ratios

    (dict) args: parser
    """
    nbins_eta, nbins_pt = len(bins['eta'])-1, len(bins['pt'])-1
    
    corr_dir = "correction_files/"
    
    mz_mc = np.loadtxt(corr_dir+f"mc{corr}_response.txt")
    res_mc = np.loadtxt(corr_dir+f"mc{corr}_resolution.txt")
    mz_dt = np.loadtxt(corr_dir+f"dt{corr}_response.txt")
    res_dt = np.loadtxt(corr_dir+f"dt{corr}_resolution.txt")
    
    ratio_mu = (91.1876+mz_dt)/(91.1876+mz_mc)
    # ratio_sigma = res_dt/res_mc * ratio_mu
    ratio_sigma = res_dt/res_mc
    labels1, labels2 = np.zeros(nbins_eta), np.zeros(nbins_pt)

    for i in range(nbins_eta):
        labels1[i] = .5*(bins['eta'][i] + bins['eta'][i+1])

    for j in range(nbins_pt):
        labels2[j] = .5*(bins['pt'][j] + bins['pt'][j+1])

    mean_min, mean_max = 0.998, 1.002
    std_min, std_max = 0.8, 1.2
    utils.plot2d(
        matrix=np.flip(ratio_mu, 0), 
        outfile=f'plots/unsmeared_massfits{corr}/mu',
        title=r'$\frac{M_Z (\mathrm{data})}{M_Z (\mathrm{mc})}$',
        x=r'$p_\mathrm{T}$',
        xbins=bins['pt'],
        y=r'$\eta$',
        ybins=bins['eta'],
        cmin=mean_min,
        cmax=mean_max,
        xticks=np.around(labels2, 1),
        yticks=labels1,
    )

    utils.plot2d(
        matrix=np.flip(ratio_sigma, 0),
        outfile=f'plots/unsmeared_massfits{corr}/sigma',
        title=r'$\frac{\mathrm{res}(M_Z(\mathrm{data}))}{\mathrm{res}(M_Z (\mathrm{mc}))}$',
        x=r'$p_\mathrm{T}$',
        xbins=bins['pt'],
        y=r'$\eta$',
        ybins=bins['eta'],
        cmin=std_min,
        cmax=std_max,
        xticks=np.around(labels2, 1),
        yticks=labels1,
    )


ROOT.gInterpreter.Declare("float gaus(){return gRandom->Gaus(0,1);}")
        

def smear_pt(datasets, ch, bins, args, corr=''):
    """
    function to smear pt w/ additional smearing factor
    
    (dict) datasets: dictionary with data samples
    (str) ch: current channel
    (dict) bins: dictionary with bin edges
    (parser) args: parser
    (str) corr: '' if not corrected and '_corr' if corrected
    """
    # prepare dataframe
    tchain = load_tchain(datasets, ch)        
    rdf = ROOT.RDataFrame(tchain)
    rdf = select_events(rdf, ch, args)

    # prepare saving
    corr_dir = "correction_files/"

    # load fit results for correction
    mz_mc = np.loadtxt(corr_dir+f"mc{corr}_response.txt")
    mz_dt = np.loadtxt(corr_dir+f"dt{corr}_response.txt")
    pt_sf = (91.1876+mz_mc) / (91.1876+mz_dt)
    res_dt = np.loadtxt(corr_dir+f"dt{corr}_resolution.txt")*pt_sf
    res_mc = np.loadtxt(corr_dir+f"mc{corr}_resolution.txt")

    # definition of smearing and scaling factors
    rdf = rdf.Define('smear_1', '0')
    rdf = rdf.Define('smear_2', '0')
    rdf = rdf.Define('scale_1', '1')
    rdf = rdf.Define('scale_2', '1')

    for e in range(len(bins['eta'])-1):
        for p in range(len(bins['pt'])-1):

            # mc resolution should not change if already larger than in data
            if res_mc[e][p] > res_dt[e][p]:
                res_sf = 1
            else:
                res_sf = res_dt[e][p] / res_mc[e][p]

            # filters for correct pt / eta bin
            f1 = f"(pt_1 > {bins['pt'][p]} && pt_1 < {bins['pt'][p+1]}) && "\
                    f"(eta_1 > {bins['eta'][e]} && eta_1 < {bins['eta'][e+1]})"
            f2 = f"(pt_2 > {bins['pt'][p]} && pt_2 < {bins['pt'][p+1]}) && "\
                    f"(eta_2 > {bins['eta'][e]} && eta_2 < {bins['eta'][e+1]})"

            # do not scale momentum in mc
            if 'mc' in ch:
                scl1 = '1'
                scl2 = '1'

            # scale momentum in data
            else:
                scl1 = 'double rsp_sf;'\
                    f'if ({f1}) rsp_sf = {pt_sf[e][p]};'\
                    f'else rsp_sf = scale_1;'\
                    'return rsp_sf;'
            
                scl2 = 'double rsp_sf;'\
                    f'if ({f2}) rsp_sf = {pt_sf[e][p]};'\
                    f'else rsp_sf = scale_2;'\
                    'return rsp_sf;'

            # smear momentum in data and mc for evaluation of smearing factor
            res1 = 'double res_sf;'\
                f'if ({f1}) res_sf = {res_mc[e][p]}/91.1876*sqrt({res_sf**2}-1)*(float)(gaus());'\
                f'else res_sf = smear_1;'\
                'return res_sf;'

            res2 = 'double res_sf;'\
                f'if ({f2}) res_sf = {res_mc[e][p]}/91.1876*sqrt({res_sf**2}-1)*(float)(gaus());'\
                f'else res_sf = smear_2;'\
                'return res_sf;'
            
            # apply calculations to calculate smearing / scaling factors
            rdf = rdf.Redefine('scale_1', scl1)
            rdf = rdf.Redefine('scale_2', scl2)
            rdf = rdf.Redefine('smear_1', res1)
            rdf = rdf.Redefine('smear_2', res2)

    # define basic output values
    outputs = ['pt_1', 'pt_2', 'eta_1', 'eta_2', 'weight']

    # apply different additional smearings, calculate outputs and save ntuples
    for i in np.linspace(0, 4, 21).round(1):
        istr = str(i).replace(".", "")
        rdf = rdf.Define("pt_1_c"+istr, f'pt_1*scale_1*(1+{i}*smear_1)')
        rdf = rdf.Define("pt_2_c"+istr, f'pt_2*scale_2*(1+{i}*smear_2)')

        rdf = calc_m(rdf, '_c'+istr)
        outputs += ["pt_1_c"+istr, "pt_2_c"+istr, "m_vis_c"+istr]

    rdf.Snapshot('ntuple', f'root_files/{ch}_smeared_ntuples.root', outputs)



def smeared_hists(ch, bins):
    """
    function that creates smeared histograms

    (str) ch: current channel
    (dict) bins: dictionary with bin edges
    (str) corr: '' if uncorrected, '_corr' if corrected
    """
    # open ntuple with smeared quantities in rdf
    tf = ROOT.TFile(f'root_files/{ch}_smeared_hists.root', 'RECREATE')
    rdf = ROOT.RDataFrame('ntuple', f'root_files/{ch}_smeared_ntuples.root')
    
    # create histograms in bins of eta and pt with different extra smearing i
    for i in tqdm(np.linspace(0,4,21).round(1)):
        istr = str(i).replace(".", "")
    
        for e in range(len(bins['eta'])-1):
            for p in range(len(bins['pt'])-1):
                
                f1 = f"((pt_1_c{istr} > {bins['pt'][p]}) &&"\
                        f"(pt_1_c{istr} < {bins['pt'][p+1]})) && "\
                        f"((eta_1 > {bins['eta'][e]}) &&"\
                        f"(eta_1 < {bins['eta'][e+1]}))"
                f2 = f1.replace('pt_1', 'pt_2').replace('eta_1', 'eta_2')

                rdf_help = rdf.Filter(f'({f1}) || ({f2})')
                    
                hist = rdf_help.Histo1D(
                    (f'm_vis_c{istr}_{e}{p}_{ch}', 'visible mass', 200, 50, 130),
                    'm_vis_c'+istr,
                    'weight'
                )
                hist.Write()
    tf.Close()
        

def fit_smeared_hists(ch, e, p):
    """
    function to fit smeared hists and plot results

    (str) ch: current channel
    (int) e: eta bin
    (int) p: pt bin
    """
    # load root file with histograms
    tf = ROOT.TFile(f'root_files/{ch}_smeared_hists.root', 'READ')

    plotdir = f"plots/smeared_massfits/{ch}"
    os.makedirs(plotdir, exist_ok=True)

    # perform fits and save resolutions and their unc.
    xyerr = []
    for i in np.linspace(0, 4, 21).round(1):
        istr = str(i).replace('.', '')
        hist = tf.Get(f'm_vis_c{istr}_{e}{p}_{ch}')

        fit_results, fit_errors, _e, _p = utils.roofit_mass(
            hist_mc = hist,
            hist_dt = False, # for technical reasons
            plot=f'{plotdir}/{e}{p}_{istr}',
        )

        xyerr.append([i, e, p, fit_results['mc'][1], fit_errors['mc'][1]])

    np.savetxt(f'fit_results/xyerr_{ch}_{e}{p}.txt', xyerr)


def plot_smearing(n_eta, n_pt):
    """
    function to plot results of fit to resolutions

    (int) n_eta: number of bins in eta
    (int) n_pt: number of bins in pt
    """

    corr_dir = "correction_files/"
    mz_mc = np.loadtxt(corr_dir+f"mc_response.txt")
    mz_dt = np.loadtxt(corr_dir+f"dt_response.txt")
    pt_sf = (91.1876+mz_mc) / (91.1876+mz_dt)
    res_dt = np.loadtxt(corr_dir+f"dt_resolution.txt")*pt_sf
    res_mc = np.loadtxt(corr_dir+f"mc_resolution.txt")

    res_sf = np.zeros((n_eta, n_pt))
    res_mc, res_dt = np.zeros((n_eta, n_pt)), np.zeros((n_eta, n_pt))

    for e in range(n_eta):
        for p in range(n_pt):
            mc = np.loadtxt(f'fit_results/xyerr_mcmm_{e}{p}.txt')
            dt = np.loadtxt(f'fit_results/xyerr_2022Cmm_{e}{p}.txt')
            #print(dt[:,4])
            x_sf, y_mc, y_dt = utils.plot_resol(
                mc[:,0],
                mc[:,3],
                mc[:,4],
                res_mc[e][p],
                dt[:,3],
                dt[:,4],
                res_dt[e][p],
                f'{e}{p}'
            )
            res_sf[e][p] = x_sf
            res_mc[e][p] = y_mc
            res_dt[e][p] = y_dt
    np.savetxt('correction_files/res_sf.txt', res_sf)
    np.savetxt('correction_files/mc_resolution_corr.txt', res_mc)
    np.savetxt('correction_files/dt_resolution_corr.txt', res_dt)
    
    
def generate_files(arguments, nthreads):
    pool = Pool(nthreads, initargs=(RLock(),), initializer=tqdm.set_lock)
    for _ in tqdm(
        pool.imap_unordered(job_wrapper_smeared_fits, arguments),
        total=len(arguments),
        desc="Total progess",
        dynamic_ncols=True,
        leave=True,
        ):
        pass


def job_wrapper_fits(args):
    return utils.roofit_mass(*args)


def job_wrapper_smeared_fits(args):
    return fit_smeared_hists(*args)


if __name__=='__main__':
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.ROOT.EnableImplicitMT(16)
    args = parse_args()
    
    basedir = '/ceph/jdriesch/CROWN_samples/Run3V07/ntuples/2022/'
    muon = 'Muon_Run2022C-PromptNanoAODv10'
    dy = 'DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv11-126X'
    datasets = {
        '2022Cmm': [
            f'{basedir}/Single{muon}/mm/Single{muon}_*.root',
            f'{basedir}/{muon}/mm/{muon}*.root',
        ],
        'mcmm': [
            f'{basedir}/{dy}/mm/{dy}*.root'
        ]
    }

    if args.corr:
        corr = '_corr'
    else:
        corr = ''
        
    # smearing_factor(datasets, bins, args)
    bins = {
        'eta': [-2.4, -1.2, 0., 1.2, 2.4],
        'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]
    }
    n_eta, n_pt = len(bins['eta'])-1, len(bins['pt'])-1

    for ch in datasets.keys():
        if args.histogram:
            make_hists(datasets, ch, bins, args, corr)

        if args.smear:
            smear_pt(datasets, ch, bins, args)
            smeared_hists(ch, bins)

    if args.fit:
        plot_fits(n_eta, n_pt, corr)
        make_ratio(n_eta, n_pt, corr)
        make2d(bins, corr)

    for ch in datasets.keys():
        if args.calculation:
            arguments = []
            nthreads = 16
            for e in range(len(bins['eta'])-1):
                for p in range(len(bins['pt'])-1):
                    arguments.append((ch, e, p))

            generate_files(arguments, nthreads)

    if args.plot:
        plot_smearing(n_eta, n_pt)
