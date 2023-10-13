import ROOT
from tqdm import tqdm
import numpy as np
import utils
import matplotlib.pyplot as plt
import scipy.optimize as opt

ROOT.gInterpreter.Declare("""
        float gaus(){
            return gRandom->Gaus(0,1);
        }
        """)


def correct(mcdt, bins = None):
    output = ['weight', 'pt_1', 'pt_2', 'eta_1', 'eta_2', 'q_1', 'q_2', 
    'trg_single_mu24_1', 'trg_single_mu24_2', 'smear_1', 'smear_2', 'run']
    
    
    basepath = '/work/jdriesch/phd/Z_early_Run3/'\
               'corrections/lepton/correction_files/Run3/mm/'
    mz_mc = np.loadtxt(basepath+'mz_mc_eta_pt_02.txt')
    mz_dt = np.loadtxt(basepath+'mz_dt_eta_pt_02.txt')
    pt_sf = (91.1876+mz_mc) / (91.1876+mz_dt)
    # mz_res_mc = np.loadtxt(basepath+'resmz_mc_eta_pt_02.txt')
    # mz_res_dt = np.loadtxt(basepath+'resmz_dt_eta_pt_02.txt') * pt_sf
    mz_res_mc = np.loadtxt('results/res_mc_corr.txt')
    mz_res_dt = np.loadtxt('results/res_dt_corr.txt') * pt_sf
    

    rdf = ROOT.RDataFrame('ntuple', mcdt['mc'])
    #rdf = ROOT.RDataFrame('ntuple', mcdt['dt'])

    weight = 'genweight*sumwWeight*crossSectionPerEventWeight*sf_trk*sf_sta*sf_id*sf_iso*sf_trg'
   
    rdf = rdf.Define('weight', weight)
    #rdf = rdf.Define('weight', '1')
    rdf = rdf.Define('smear_1', '0')
    rdf = rdf.Define('smear_2', '0')

    if bins:
        filter = '(pt_{n} > {p_l} && pt_{n} < {p_r}) \
                   && (eta_{n} > {e_l} && eta_{n} < {e_r})'

        for p in range(len(bins['pt'])-1):
            for e in range(len(bins['eta'])-1):

                if mz_res_mc[e][p] > mz_res_dt[e][p]:
                    res_sf = 1
                else:
                    res_sf = mz_res_dt[e][p]/mz_res_mc[e][p]

                def_temp = 'double res_sf;'\
                           'if ({f}) res_sf = {res}/91.1876*sqrt({sf}-1)*(float)(gaus());'\
                           'else res_sf = smear_{n};'\
                           'return res_sf;'

                rdf = rdf.Redefine(
                    'smear_1', 
                    def_temp.format(
                        f = filter.format(
                            n=1, 
                            p_l=bins['pt'][p], p_r=bins['pt'][p+1], 
                            e_l=bins['eta'][e], e_r=bins['eta'][e+1]
                            ), 
                        res=mz_res_mc[e][p],
                        sf=res_sf**2, 
                        n=1
                        )
                    )
                rdf = rdf.Redefine(
                    'smear_2', 
                    def_temp.format(
                        f = filter.format(
                            n=2, 
                            p_l=bins['pt'][p], p_r=bins['pt'][p+1], 
                            e_l=bins['eta'][e], e_r=bins['eta'][e+1]
                            ), 
                        res=mz_res_mc[e][p],
                        sf=res_sf**2, 
                        n=2
                        )
                    )
        print('finished calculating smearing factors')

    else:
        rdf = rdf.Redefine('smear_1', '1./pt_1 * (float)(gaus())')
        rdf = rdf.Redefine('smear_2', '1./pt_2 * (float)(gaus())')
    
    for i in np.linspace(0,2,11).round(1):
        istr = str(i).replace(".", "")
        rdf = rdf.Define('pt_1_c'+str(istr), 'pt_1*(1+{}*smear_1)'.format(i))
        rdf = rdf.Define('pt_2_c'+str(istr), 'pt_2*(1+{}*smear_2)'.format(i))

        rdf = rdf.Define(
            'lv_1'+str(istr), 
            'ROOT::Math::PtEtaPhiMVector p(pt_1_c{}, eta_1, phi_1, mass_1); \
                return p'.format(istr)
            )
        rdf = rdf.Define(
            'lv_2'+str(istr), 
            'ROOT::Math::PtEtaPhiMVector p(pt_2_c{}, eta_2, phi_2, mass_2); \
                return p'.format(istr)
            )
        
        rdf = rdf.Define('dimuon'+str(istr), 'lv_1{i} + lv_2{i}'.format(i = istr))
        rdf = rdf.Define('pt_vis_c'+str(istr), 'dimuon{}.Pt()'.format(istr))
        rdf = rdf.Define('m_vis_c'+str(istr), 'dimuon{}.M()'.format(istr))

        output.append('m_vis_c'+str(istr))
    print('finished applying smearing factors')
    #'''
    rdf.Snapshot(
        'ntuple', 
        '/ceph/jdriesch/CROWN_samples/Run3V02/resol_impact/mc_corr_fac_sf_it.root',
        output
        )    
    '''
    print(output)
    rdf.Snapshot(
        'ntuple', 
        '/ceph/jdriesch/CROWN_samples/Run3V02/resol_impact/dt_corr_fac_sf.root',
        output
        )
    '''
    print('corrections applied')


def make_hists(mcdt, bins):
    nbins_pt, nbins_eta = len(bins['pt']), len(bins['eta'])

    rdf = ROOT.RDataFrame('ntuple', mcdt['mc'])
    #rdf = ROOT.RDataFrame('ntuple', mcdt['dt'])

    event_selection = '(q_1*q_2 < 0) && (trg_single_mu24_1 || trg_single_mu24_2)'
    event_selection += f' && (run >= 355862 || run == 1)'
    event_selection += f' && (run <= 357482 || run == 1)'

    rdf = rdf.Filter(event_selection)
    #rdf_dt = rdf_dt.Filter(event_selection)    

    filter = '((pt_{n} > {p_l} && pt_{n} < {p_r}) \
                && (eta_{n} > {e_l} && eta_{n} < {e_r}))'

    for p in tqdm(range(nbins_pt-1)):
        for e in range(nbins_eta-1):
            f1 = filter.format(
                n=1, 
                p_l=bins['pt'][p], p_r=bins['pt'][p+1], 
                e_l=bins['eta'][e], e_r=bins['eta'][e+1]
                )
            f2 = filter.format(
                n=2, 
                p_l=bins['pt'][p], p_r=bins['pt'][p+1], 
                e_l=bins['eta'][e], e_r=bins['eta'][e+1]
                )
            rdf_help = rdf.Filter('{} || {}'.format(f1, f2))

            file0 = ROOT.TFile('hists/mc_it/hists_corr_{}_{}_sf.root'.format(p,e), 'RECREATE')
            for i in np.linspace(0,2,11).round(1):
                istr = str(i).replace('.', '')
                hist = rdf_help.Histo1D(
                    ('m_vis_c'+str(istr), 'visible mass', 200, 50, 130), 
                    'm_vis_c'+str(istr), 
                    'weight'
                    )
                hist.Write()
            file0.Close()


def fit_hists(p,e):
    ROOT.RooMsgService.instance().setSilentMode(True)
    file0 = ROOT.TFile('hists/mc_it/hists_corr_{}_{}_sf.root'.format(p,e))
    #file0 = ROOT.TFile('hists/dt/hists_corr_{}_{}_sf.root'.format(p,e))
    file1 = ROOT.TFile('../../hists/Run3/mm/eta_pt_02/eta_{}_pt_{}.root'.format(e,p))
    hist_dt = file1.Get('data')

    # introduce fit quantities
    x = ROOT.RooRealVar('x', 'm_vis (GeV)', 85, 97)
    #x = ROOT.RooRealVar('x', 'm_vis (GeV)', 50, 130)
    x.setBins(10000,'cache')
    x.setMin('cache',0)
    x.setMax('cache',500)
    Z_mass = ROOT.RooRealVar('Z_mass', 'Z_mass', 91.1876, 60, 120)
    Z_width = ROOT.RooRealVar('Z_width', 'Z_widthan', 2.4952, 0, 10)
    Z_mass.setConstant(True)
    Z_width.setConstant(True)

    mean = ROOT.RooRealVar('mean', 'mean', 0, -10, 10)
    #sigma = ROOT.RooRealVar('sigma', 'sigma', 2, 0, 20)
    sigma = ROOT.RooRealVar('sigma', 'sigma', 1.42, 0, 20)
    n_L = ROOT.RooRealVar('n_L', 'n_L', 5, 0, 1000)
    n_R = ROOT.RooRealVar('n_R', 'n_R', 5, 0, 1000)
    alpha_L = ROOT.RooRealVar('alpha_L', 'alpha_L', 1.3, 1, 5)
    alpha_R = ROOT.RooRealVar('alpha_R', 'alpha_R', 1.2, 1, 5)

    bw = ROOT.RooBreitWigner('bw', 'BreitWigner', x, Z_mass, Z_width)
    cb = ROOT.RooCrystalBall('cb', 'CrystalBall', x, mean, sigma,
                            sigma, alpha_L, n_L, alpha_R, n_R)

    func = ROOT.RooFFTConvPdf('func', 'func', x, bw, cb)

    smear = np.linspace(0,2,11).round(1)
    resol = []
    err = []

    for i in tqdm(smear):
        istr = str(i).replace(".", "")
        sigma.setVal(max(1,1+.3*i))

        hist_mc = file0.Get('m_vis_c'+str(istr))
        hist_mc.Scale(hist_dt.Integral()/hist_mc.Integral())
        #hist.Draw()
        roohist = ROOT.RooDataHist('mc', 'mc hist', ROOT.RooArgSet(x), hist_mc)
        fitResult = func.fitTo(
            roohist, 
            ROOT.RooFit.AsymptoticError(True), 
            ROOT.RooFit.PrintEvalErrors(-1)
            )


        c1 = ROOT.TCanvas( 'c1', 'The Fit Canvas', 200, 10, 700, 500)
        c1.SetGridx()
        c1.SetGridy()
        c1.GetFrame().SetFillColor( 21 )
        c1.GetFrame().SetBorderMode(-1 )
        c1.GetFrame().SetBorderSize( 5 )
        frame = x.frame()
        frame.SetTitle('smeared by factor {}'.format(i))

        roohist.plotOn(
            frame, 
            ROOT.RooFit.DrawOption('B'), 
            ROOT.RooFit.FillStyle(0), 
            ROOT.RooFit.FillColor(ROOT.kBlue)
            )
        func.plotOn(
            frame, 
            ROOT.RooFit.LineColor(ROOT.kBlue)
            )
        chi2 = frame.chiSquare(6)

        print('Smearing by factor of {} leads to sigma of {}'.format(i, sigma.getVal()))

        frame.Draw()
        c1.Update()
        ROOT.gStyle.SetGridColor(ROOT.kGray+1)

        cmsTex=ROOT.TLatex()
        cmsTex.SetTextFont(42)
        cmsTex.SetTextSize(0.025)
        cmsTex.SetNDC()

        cmsTex.SetTextSize(0.035)
        cmsTex.DrawLatex(0.1,0.92,'#bf{CMS} #it{Preliminary}')

        cmsTex.SetTextSize(0.025)
        cmsTex.SetLineWidth(2)
        cmsTex.SetTextFont(42)

        stats = ROOT.TPaveText(0.65, 0.65, 0.88, 0.88, 'br ARC NDC')
        stats.AddText(
            'M(Z) = {} ({})'.format(
                round(mean.getVal(),3), round(mean.getAsymErrorHi(),3)
                )
            )
        stats.GetListOfLines().Last().SetTextColor(ROOT.kBlue)
        stats.AddText(
            'res(M) = {} ({})'.format(
                round(sigma.getVal(),3), round(sigma.getAsymErrorHi(),3)
                )
            )
        stats.GetListOfLines().Last().SetTextColor(ROOT.kBlue)
        stats.AddText('chi2/dof = {}'.format(round(chi2,2)))
        stats.GetListOfLines().Last().SetTextColor(ROOT.kBlue)

        stats.Draw('SAME')
        c1.Update()
        plot = 'results/mc_it/sf_{}_{}/smearing_factor_{}'.format(
            p,e,str(i).replace(".", "")
            )
        c1.SaveAs('{}.png'.format(plot))
        c1.SaveAs('{}.pdf'.format(plot))
        
        resol.append(sigma.getVal())
        err.append(sigma.getAsymErrorHi())

    xyerr = [smear, resol, err]
    print(xyerr)

    np.savetxt('results/mc_it/sf_{}_{}/xyerr.txt'.format(p,e), xyerr)
    return

def plot_resol(p,e):
    infile = np.loadtxt('results/mc_it/sf_{}_{}/xyerr.txt'.format(p,e))
    basepath = '/work/jdriesch/phd/Z_early_Run3/corrections/lepton/correction_files/Run3/mm/'
    mz_mc = np.loadtxt(basepath+'mz_mc_eta_pt_02.txt')
    mz_dt = np.loadtxt(basepath+'mz_dt_eta_pt_02.txt')
    pt_sf = (91.1876+mz_mc) / (91.1876+mz_dt)
    #mz_res_mc = np.loadtxt(basepath+'resmz_mc_eta_pt_02.txt') 
    #mz_res_dt = np.loadtxt(basepath+'resmz_dt_eta_pt_02.txt') * pt_sf
    mz_res_mc = np.loadtxt('results/res_mc_corr.txt')
    mz_res_dt = np.loadtxt('results/res_dt_corr.txt') * pt_sf

    xdata = infile[0]
    ydata = infile[1]
    yerr  = infile[2]
    y_dt  = mz_res_dt[e][p]
    y_mc  = mz_res_mc[e][p]

    plt.errorbar(xdata, ydata, yerr, label='data', marker='+', linestyle="")
    eps=0.0000000001
    par, pcov = opt.curve_fit(
        func, 
        xdata, 
        ydata, 
        bounds=((0, min(ydata[0], y_mc)), (np.inf, max(ydata[0], y_mc)))
        )
    # par, pcov = opt.curve_fit(func, xdata, ydata, bounds=((0, min(ydata[0], y_dt)), (np.inf, max(ydata[0], y_dt))))

    x = np.linspace(0, 2.4, 1000)
    y = func(x, *par)
    if y[-1] < y_dt:
        x_sf = 1
    else:
        x_sf = x[np.where(y > y_dt)[0][0]]
    print(x_sf, y[0])
    plt.plot(x, y, label='fit')
    plt.plot(x, np.ones_like(x)*y_dt, label='data resolution')
    plt.plot(x, np.ones_like(x)*y_mc, label='mc resolution')
    plt.plot(x_sf, y_dt, label='additional smearing', marker='x', linestyle="")
    plt.legend()

    plt.text(0.1, 3, 'smearing: {}'.format(round(x_sf,3)))
    plt.text(0.1, 2.7, 'gradient: {}'.format(round(par[0],3)))
    plt.xlabel('smearing factor of pt')
    plt.ylabel('additional smearing of m_vis')
    plt.savefig('results/mc_it/fits/resol_{}_{}.png'.format(p,e))
    plt.clf()
    return x_sf, y[0]


def func(x, m, x0):
     return np.sqrt((m * x)**2 + x0**2)


if __name__=='__main__':
    mcdt = {
        'mc': '/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC/2022/'\
              'DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/'\
              'DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root',
        'dt': [
            '/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC/2022/'\
              'Muon_Run2022C-PromptReco-v1/mm/Muon_Run2022C-PromptReco-v1_*.root',
            '/ceph/moh/CROWN_samples/Run3V03/ntuples_xsec_sf_EraC/2022/'\
              'SingleMuon_Run2022C-PromptReco-v1/mm/SingleMuon_Run2022C-PromptReco-v1_*.root',
        ]
    }
    bins = {
        'eta': [-2.4, -1.2, 0., 1.2, 2.4], 
        'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]
        }
    #correct(mcdt, bins)

    mcdt_2 = {
        'mc': '/ceph/jdriesch/CROWN_samples/Run3V02/resol_impact/mc_corr_fac_sf_it.root',
        'dt': '/ceph/jdriesch/CROWN_samples/Run3V02/resol_impact/dt_corr_fac_sf.root',
    }
    #make_hists(mcdt_2, bins)

    fit = False
    plot = True
    sf = np.zeros((len(bins['eta'])-1, len(bins['pt'])-1))
    res = np.zeros((len(bins['eta'])-1, len(bins['pt'])-1))
    for p in range(len(bins['pt'])-1):
        for e in range(len(bins['eta'])-1):
            utils.usedir('results/mc_it/sf_{}_{}/'.format(p,e), overwrite=True)
            if fit:# and p>2 and e==1:
                fit_hists(p,e)
            #if plot and p==3 and e==0:
            sf[e][p], res[e][p] = plot_resol(p,e)
    if plot:
        np.savetxt('results/sf_extra_mc_it.txt', sf)
        np.savetxt('results/res_mc_corr_it.txt', res)
        print("")