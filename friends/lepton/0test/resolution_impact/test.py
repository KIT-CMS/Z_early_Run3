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

def make_test_file(mcdt, bins):
    output = ['weight', 'pt_1', 'pt_2', 'eta_1', 'eta_2', 'q_1', 'q_2', 'trg_single_mu24_1', 'trg_single_mu24_2', 'smear_1', 'smear_2']

    basepath = '/work/jdriesch/phd/Z_early_Run3/'\
               'corrections/lepton/correction_files/Run3/mm/'
    mz_mc = np.loadtxt(basepath+'mz_mc_eta_pt_02.txt')
    mz_dt = np.loadtxt(basepath+'mz_dt_eta_pt_02.txt')
    pt_sf = (91.1876+mz_mc) / (91.1876+mz_dt)
    mz_res_mc = np.loadtxt(basepath+'resmz_mc_eta_pt_02.txt')
    mz_res_dt = np.loadtxt(basepath+'resmz_dt_eta_pt_02.txt') * pt_sf

    rdf = ROOT.RDataFrame('ntuple', mcdt['mc'])

    event_selection = '(q_1*q_2 < 0) && (trg_single_mu24_1 || trg_single_mu24_2)'
    rdf = rdf.Filter(event_selection)

    filter = '((pt_{n} > {p_l} && pt_{n} < {p_r}) \
                && (eta_{n} > {e_l} && eta_{n} < {e_r}))'

    p, e = 0, 0
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
    rdf = rdf.Filter('{} || {}'.format(f1, f2))



    weight = 'genweight*sumwWeight*crossSectionPerEventWeight*sf_trk*sf_sta*sf_id*sf_iso*sf_trg'
   
    rdf = rdf.Define('weight', weight)
    rdf = rdf.Define('smear_1', '0')
    rdf = rdf.Define('smear_2', '0')


    for p in range(len(bins['pt'])-1):
        for e in range(len(bins['eta'])-1):

            if mz_res_mc[e][p] > mz_res_dt[e][p]:
                res_sf = 1
                print('res_data < res_mc for bin [pt, eta] = [{},{}]'.format(p,e))
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

    
    for i in np.linspace(0,2,11):
        istr = str(int(i*10))
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
        output.append('pt_1_c'+str(istr))
        output.append('pt_2_c'+str(istr))
        output.append('pt_vis_c'+str(istr))

    print('finished applying smearing factors')

    rdf.Snapshot(
        'ntuple', 
        'mc_sf.root',
        output
        )


def test():
    rdf = ROOT.RDataFrame('ntuple', 'mc_sf.root')
    h1 = rdf.Histo2D(('h1', 'm_vis_c6 vs smear_1', 100, 60,120, 100, -2, 2), 'm_vis_c6', 'smear_2')

    h2 = rdf.Histo1D(('h2', 'm_vis_c6', 100, 60, 120), 'm_vis_c6')
    h1.SetMaximum(100)
    c = ROOT.TCanvas('c', 'smearing of pt', 800, 1200)
    c.Divide(1,2)
    c.cd(1)
    h1.Draw("COLZ")
    c.cd(2)
    h2.Draw()
    #print(h1.GetMaximumBin(), h1.GetMaximum())
    c.SaveAs('smearing.pdf')

if __name__=='__main__':
    mcdt = {
        'mc': '/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf/2022/'\
              'DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/'\
              'DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root',
        'dt': '/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf_EraC_lep_corr_01_x0p60/2022/'\
              'Muon_Run2022C-PromptReco-v1/mm/Muon_Run2022C-PromptReco-v1_*.root'
    }
    bins = {
        'eta': [-2.4, -1.2, 0., 1.2, 2.4], 
        'pt': [25, 30, 35, 39, 42, 45, 49, 60, 120, 9999]
        }
    make_test_file(mcdt, bins)
    test()
