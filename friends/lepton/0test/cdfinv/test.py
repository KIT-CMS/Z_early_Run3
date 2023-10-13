import ROOT
import argparse
import os
from tqdm import tqdm
import numpy as np
import array as a
import json
import glob
import utils

def makehists(mcdt, corr=''):
    df_mc = ROOT.RDataFrame('ntuple', mcdt['mc'])
    df_dt = ROOT.RDataFrame('ntuple', mcdt['dt'])

    weight = "genweight*sumwWeight*crossSectionPerEventWeight*sf_sel*sf_trg"

    df_mc = df_mc.Define("weight", weight)
    file0 = ROOT.TFile('hists{}.root'.format(corr), 'RECREATE')

    eta = [-2.4, -1.2, 0, 1.2, 2.4]
    for i in range(len(eta)-1):
        pt_1_mc = df_mc.Filter('eta_1 > '+str(eta[i])).Filter('eta_1 < '+str(eta[i+1])).Histo1D(('pt_1_mc_eta_'+str(i), 'pt_1', 1000, 20, 250), 'pt_1', 'weight')
        pt_1_dt = df_dt.Filter('eta_1 > '+str(eta[i])).Filter('eta_1 < '+str(eta[i+1])).Histo1D(('pt_1_dt_eta_'+str(i), 'pt_1', 1000, 20, 250), 'pt_1')
        pt_2_mc = df_mc.Filter('eta_2 > '+str(eta[i])).Filter('eta_2 < '+str(eta[i+1])).Histo1D(('pt_2_mc_eta_'+str(i), 'pt_2', 1000, 20, 250), 'pt_2', 'weight')
        pt_2_dt = df_dt.Filter('eta_2 > '+str(eta[i])).Filter('eta_2 < '+str(eta[i+1])).Histo1D(('pt_2_dt_eta_'+str(i), 'pt_2', 1000, 20, 250), 'pt_2')

        pt_1_mc.Write()
        pt_1_dt.Write()        
        pt_2_mc.Write()
        pt_2_dt.Write()

    file0.Close()


def plotmass(mcdt, corr=''):
    df_mc = ROOT.RDataFrame('ntuple', mcdt['mc'])
    df_dt = ROOT.RDataFrame('ntuple', mcdt['dt'])    
    weight = "genweight*sumwWeight*crossSectionPerEventWeight*sf_sel*sf_trg"

    df_mc = df_mc.Define("weight", weight)

    m_vis_mc = df_mc.Histo1D(('m_vis_mc', 'm_vis', 100, 60, 120), 'm_vis', 'weight')
    m_vis_dt = df_dt.Histo1D(('m_vis_dt', 'm_vis', 100, 60, 120), 'm_vis')    
    m_vis_corr_mc = df_mc.Histo1D(('m_vis_corr_mc', 'm_vis', 100, 60, 120), 'm_vis_corr', 'weight')
    m_vis_corr_dt = df_dt.Histo1D(('m_vis_corr_dt', 'm_vis', 100, 60, 120), 'm_vis_corr')
    m_vis_c_mc = df_mc.Histo1D(('m_vis_c_mc', 'm_vis', 100, 60, 120), 'm_vis_c', 'weight')

    file0 = ROOT.TFile('m_hists.root', 'RECREATE')
    m_vis_mc.Write()
    m_vis_dt.Write()
    m_vis_corr_mc.Write()
    m_vis_corr_dt.Write()
    m_vis_c_mc.Write()
    file0.Close()


    file0 = ROOT.TFile('m_hists.root')
    m_vis_mc = file0.Get('m_vis_mc')
    m_vis_dt = file0.Get('m_vis_dt')
    m_vis_corr_mc = file0.Get('m_vis_corr_mc')
    m_vis_corr_dt = file0.Get('m_vis_corr_dt')
    m_vis_c_mc = file0.Get('m_vis_c_mc')

    m_vis_mc.Scale(m_vis_dt.Integral()/m_vis_mc.Integral())

    utils.plot_ratio(plots = {'mc': m_vis_mc, 'dt': m_vis_dt}, rcolors = {'mc': ROOT.kBlue, 'dt': ROOT.kBlack}, title='cdf inversion method', outfile='plots/m_vis.pdf', text=['','',''], evts=-1, xrange=[80,102])
    """
    c = ROOT.TCanvas('c', 'c', 1800, 600)
    c.Divide(3,1)
    c.cd(1)
    m_vis_mc.Draw()
    m_vis_dt.Draw('same')
    c.cd(2)
    m_vis_corr_mc.Draw()
    m_vis_corr_dt.Draw('same')
    c.cd(3)
    m_vis_c_mc.Draw()
    m_vis_dt.Draw('same')
    c.SaveAs('plots/m_vis_comp.pdf')
    c.SaveAs('plots/m_vis_comp.png')
    """



def printhists(corr=''):
    file0 = ROOT.TFile('hists{}.root'.format(corr))
    eta = [-2.4, -1.2, 0, 1.2, 2.4]
    for i in range(len(eta)-1):
        pt_1_mc = file0.Get('pt_1_mc_eta_'+str(i))
        pt_1_dt = file0.Get('pt_1_dt_eta_'+str(i))        
        pt_2_mc = file0.Get('pt_2_mc_eta_'+str(i))
        pt_2_dt = file0.Get('pt_2_dt_eta_'+str(i))
        #print(pt_1_mc.Integral())

        pt_1_mc.Scale(pt_1_dt.Integral()/pt_1_mc.Integral())
        pt_2_mc.Scale(pt_2_dt.Integral()/pt_2_mc.Integral())

        c = ROOT.TCanvas('c', 'c', 800, 600)
        pt_1_mc.Draw()
        pt_1_dt.Draw('same')
        c.SaveAs('plots/pt_1_eta_{}{}.pdf'.format(i, corr))
        c.SaveAs('plots/pt_1_eta_{}{}.png'.format(i, corr))        
        
        c = ROOT.TCanvas('c', 'c', 800, 600)
        pt_2_mc.Draw()
        pt_2_dt.Draw('same')
        c.SaveAs('plots/pt_2_eta_{}{}.pdf'.format(i, corr))
        c.SaveAs('plots/pt_2_eta_{}{}.png'.format(i, corr))


def makepdfs(corr=''):
    file0 = ROOT.TFile('hists{}.root'.format(corr))
    eta = [-2.4, -1.2, 0, 1.2, 2.4]
    ws = ROOT.RooWorkspace('ws')

    for i in range(len(eta)-1):
        pt_1_mc = file0.Get('pt_1_mc_eta_'+str(i))
        pt_1_dt = file0.Get('pt_1_dt_eta_'+str(i))        
        pt_2_mc = file0.Get('pt_2_mc_eta_'+str(i))
        pt_2_dt = file0.Get('pt_2_dt_eta_'+str(i))

        x = ROOT.RooRealVar("x", "m_vis (GeV)", 20, 130)
        x.setBins(10000,"cache")
        x.setMin("cache",0)
        x.setMax("cache",500)


        pt_1_mc_h = ROOT.RooDataHist('pt_1_mc_eta_{}_h'.format(i), 'mc', ROOT.RooArgSet(x), pt_1_mc)
        pt_2_mc_h = ROOT.RooDataHist('pt_2_mc_eta_{}_h'.format(i), 'mc', ROOT.RooArgSet(x), pt_2_mc)
        pt_1_mc_pdf  = ROOT.RooHistPdf('pt_1_mc_eta_{}_pdf'.format(i), 'mc', ROOT.RooArgSet(x), pt_1_mc_h)
        pt_2_mc_pdf  = ROOT.RooHistPdf('pt_2_mc_eta_{}_pdf'.format(i), 'mc', ROOT.RooArgSet(x), pt_2_mc_h)

        pt_1_dt_h = ROOT.RooDataHist('pt_1_dt_eta_{}_h'.format(i), 'dt', ROOT.RooArgSet(x), pt_1_dt)
        pt_2_dt_h = ROOT.RooDataHist('pt_2_dt_eta_{}_h'.format(i), 'dt', ROOT.RooArgSet(x), pt_2_dt)
        pt_1_dt_pdf  = ROOT.RooHistPdf('pt_1_dt_eta_{}_pdf'.format(i), 'dt', ROOT.RooArgSet(x), pt_1_dt_h)
        pt_2_dt_pdf  = ROOT.RooHistPdf('pt_2_dt_eta_{}_pdf'.format(i), 'dt', ROOT.RooArgSet(x), pt_2_dt_h)

        ws.Import(pt_1_mc_pdf)
        ws.Import(pt_2_mc_pdf)
        ws.Import(pt_1_dt_pdf)
        ws.Import(pt_2_dt_pdf)
    ws.writeToFile('ws{}.root'.format(corr))


def plotcdfs(corr=''):
    ws = ROOT.TFile('ws{}.root'.format(corr)).Get('ws')
    x = ws.var('x')
    eta = [-2.4, -1.2, 0, 1.2, 2.4]

    for i in range(len(eta)-1):
        cdf_mc1 = ws.pdf('pt_1_mc_eta_{}_pdf'.format(i)).createCdf(x)
        cdf_mc2 = ws.pdf('pt_2_mc_eta_{}_pdf'.format(i)).createCdf(x)
        cdf_dt1 = ws.pdf('pt_1_dt_eta_{}_pdf'.format(i)).createCdf(x)
        cdf_dt2 = ws.pdf('pt_2_dt_eta_{}_pdf'.format(i)).createCdf(x)    

        frame = x.frame()
        frame.SetTitle('')
        cdf_mc1.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGreen))
        cdf_dt1.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlack))

        c1 = ROOT.TCanvas( 'c1', 'The Fit Canvas', 200, 10, 700, 500 )
        c1.SetGridx()
        c1.SetGridy()
        c1.GetFrame().SetFillColor( 21 )
        c1.GetFrame().SetBorderMode(-1 )
        c1.GetFrame().SetBorderSize( 5 )

        frame.Draw()
        c1.Update()
        ROOT.gStyle.SetGridColor(ROOT.kGray+1)

        c1.Update()
        name = 'cdf_1_eta_{}{}'.format(i, corr)
        c1.SaveAs("plots/{}.png".format(name))
        c1.SaveAs("plots/{}.pdf".format(name))
        

        frame = x.frame()
        frame.SetTitle('')
        cdf_mc2.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGreen))
        cdf_dt2.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlack))
        c2 = ROOT.TCanvas( 'c2', 'The Fit Canvas', 200, 10, 700, 500 )
        c2.SetGridx()
        c2.SetGridy()
        c2.GetFrame().SetFillColor( 21 )
        c2.GetFrame().SetBorderMode(-1 )
        c2.GetFrame().SetBorderSize( 5 )

        frame.Draw()
        c2.Update()
        ROOT.gStyle.SetGridColor(ROOT.kGray+1)

        c2.Update()
        name = 'cdf_2_eta_{}{}'.format(i, corr)
        c2.SaveAs("plots/{}.png".format(name))
        c2.SaveAs("plots/{}.pdf".format(name))


def correct(mcdt):
    eta = [-2.4, -1.2, 0, 1.2, 2.4]
    cdf = {}

    ws = ROOT.TFile('ws.root').Get('ws')
    ws.Print()
    x = ws.var('x')

    for i in range(len(eta)-1):
        cdf['pt_1_mc_eta_'+str(i)] = ws.pdf('pt_1_mc_eta_{}_pdf'.format(i)).createCdf(x)
        cdf['pt_2_mc_eta_'+str(i)] = ws.pdf('pt_2_mc_eta_{}_pdf'.format(i)).createCdf(x)
        cdf['pt_1_dt_eta_'+str(i)] = ws.pdf('pt_1_dt_eta_{}_pdf'.format(i)).createCdf(x)
        cdf['pt_2_dt_eta_'+str(i)] = ws.pdf('pt_2_dt_eta_{}_pdf'.format(i)).createCdf(x)


    df_mc = ROOT.RDataFrame('ntuple', mcdt['mc'])
    df_mc = df_mc.Define('pt_1_c', '1').Define('pt_2_c', '1')

    np_rdf = df_mc.AsNumpy(['pt_1', 'pt_2', 'pt_1_corr', 'pt_2_corr', 'pt_1_c', 'pt_2_c', 'eta_1', 'eta_2',  'm_vis', 'm_vis_corr', 'phi_1', 'phi_2', 'mass_1', 'mass_2', 'genweight', 'sumwWeight', 'crossSectionPerEventWeight', 'sf_sel', 'sf_trg'])
    pt_1, pt_2 = np_rdf['pt_1'], np_rdf['pt_2']
    eta_1, eta_2 = np_rdf['eta_1'], np_rdf['eta_2']
    pt_1_c, pt_2_c = pt_1, pt_2

    for i in tqdm(range(len(pt_1))):
        # lepton 1
        x.setVal(pt_1[i])

        for j in range(len(eta)-1):
            if eta_1[i] > eta[j] and eta_1[i] < eta[j+1]:
                cdf_val = cdf['pt_1_mc_eta_'+str(j)].getVal() 

                inv_val = cdf['pt_1_dt_eta_'+str(j)].findRoot(x, x.getMin(), x.getMax(), cdf_val)
                pt_1_c[i] = inv_val

        # lepton 2
        x.setVal(pt_2[i])

        for j in range(len(eta)-1):
            if eta_2[i] > eta[j] and eta_2[i] < eta[j+1]:
                cdf_val = cdf['pt_2_mc_eta_'+str(j)].getVal() 

                inv_val = cdf['pt_2_dt_eta_'+str(j)].findRoot(x, x.getMin(), x.getMax(), cdf_val)
                pt_2_c[i] = inv_val      

    np_rdf['pt_1_c'] = pt_1_c
    np_rdf['pt_2_c'] = pt_2_c

    rdf = ROOT.RDF.MakeNumpyDataFrame(np_rdf)

    rdf = rdf.Define("lv_1", "ROOT::Math::PtEtaPhiMVector p(pt_1_c, eta_1, phi_1, mass_1); return p")
    rdf = rdf.Define("lv_2", "ROOT::Math::PtEtaPhiMVector p(pt_2_c, eta_2, phi_2, mass_2); return p")
    
    rdf = rdf.Define("dimuon", "lv_1 + lv_2")
    rdf = rdf.Define("pt_vis_c", "dimuon.Pt()")
    rdf = rdf.Define("m_vis_c", "dimuon.M()")

    rdf.Snapshot('ntuple', 'mc_corr.root', ['pt_1', 'pt_2', 'pt_1_corr', 'pt_2_corr', 'pt_1_c', 'pt_2_c', 'eta_1', 'eta_2', 'genweight', 'sumwWeight', 'crossSectionPerEventWeight', 'sf_sel', 'sf_trg', 'm_vis', 'm_vis_corr', 'm_vis_c'])





if __name__=='__main__':
    mcdt = {
        'mc': '/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf_EraC_lep_corr_01_x0p60/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_1*.root',
        'dt': '/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf_EraC_lep_corr_01_x0p60/2022/Muon_Run2022C-PromptReco-v1/mm/Muon_Run2022C-PromptReco-v1_1*.root'
    }
    #makehists(mcdt)
    #printhists()
    #makepdfs()
    #plotcdfs()
    mcdt_tocorr = {
        'mc': '/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf_EraC_lep_corr_01_x0p60/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*0.root',
    }
    #correct(mcdt_tocorr)

    mcdt_2 = {
        'mc': 'mc_corr.root',
        'dt': '/ceph/moh/CROWN_samples/Run3V02/ntuples_xsec_sf_EraC_lep_corr_01_x0p60/2022/Muon_Run2022C-PromptReco-v1/mm/Muon_Run2022C-PromptReco-v1_1*.root'
    }
    #makehists(mcdt_2, '_corr')
    #printhists('_corr')    
    #makehists(mcdt)
    #printhists()
    plotmass(mcdt_2)