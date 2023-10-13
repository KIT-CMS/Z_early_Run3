import ROOT
import argparse
import os
from tqdm import tqdm
import numpy as np
import utils as utils
import array as a
import json
import glob

bins = {'eta': [-2.4, -1.2, 0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 1500]}

bin1, bin2 = list(bins.keys())[0], list(bins.keys())[1]
nbins1, nbins2 = len(bins[bin1])-1, len(bins[bin2])-1


corr=""
fs='ee'
run='3'
v='corr_01_x06'
#v='v4'

path_corr = "correction_files/Run{run}/{fs}/{res}mz_{mcdt}_eta_pt_{v}{corr}.txt"
outdir = 'plots/Run{run}/{fs}/eta_pt{corr}_{v}/'.format(run=run, fs=fs, corr=corr, v=v)


"""
mz_mc, mz_dt = np.loadtxt(path_corr.format(run=run, fs=fs, res="", mcdt='mc', corr=corr, v=v)), np.loadtxt(path_corr.format(run=run, fs=fs, res="", mcdt='dt', corr=corr, v=v))
pt_sf = (91.1876+mz_mc) / (91.1876+mz_dt)
mz_res_mc, mz_res_dt = np.loadtxt(path_corr.format(run=run, fs=fs, res="res", mcdt='mc', corr=corr, v=v)), np.loadtxt(path_corr.format(run=run, fs=fs, res="res", mcdt='dt', corr=corr, v=v)) * pt_sf
res_sf = mz_res_dt/mz_res_mc
"""
"""
ratio_mu = (91.1876 + mz_dt)/(91.1876 + mz_mc)
ratio_sigma = mz_res_dt/mz_res_mc * ratio_mu
labels1, labels2 = np.zeros(nbins1), np.zeros(nbins2)

for i in range(nbins1):
    labels1[i] = .5*(bins[bin1][i] + bins[bin1][i+1])

for j in range(nbins2):
    labels2[j] = .5*(bins[bin2][j] + bins[bin2][j+1])

print("done")

utils.plot2d(ratio_mu, '{}mu{}.pdf'.format(outdir,corr), r'$\frac{M_Z (\mathrm{data})}{M_Z (\mathrm{mc})}$', r'$p_\mathrm{T}$', bins[bin2], r'$\eta$', bins[bin1], 0.98, 1.02, np.around(labels2, 1), labels1)
utils.plot2d(ratio_sigma, '{}sigma{}.pdf'.format(outdir,corr), r'$\frac{\mathrm{res}(M_Z(\mathrm{data}))}{\mathrm{res}(M_Z (\mathrm{mc}))}$', r'$p_\mathrm{T}$', bins[bin2], r'$\eta$', bins[bin1], 0.7, 1.3, np.around(labels2,1), labels1)
"""

path_mc = '/ceph/jdriesch/CROWN_samples/Run3V02/ntuples_xsec_sf_mu_{v}/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/*.root'.format(v=v)
#path_mc = '/ceph/jdriesch/CROWN_samples/Run3V02/ntuples_xsec_sf_mu_{v}/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/ee/*.root'.format(v=v)
path_dt = [
    '/ceph/jdriesch/CROWN_samples/Run3V02/ntuples_xsec_sf_mu_{v}/2022/SingleMuon_Run2022C-PromptReco-v1/mm/*.root'.format(v=v),
    '/ceph/jdriesch/CROWN_samples/Run3V02/ntuples_xsec_sf_mu_{v}/2022/Muon_Run2022C-PromptReco-v1/mm/*.root'.format(v=v),
]
#path_dt = ['/ceph/jdriesch/CROWN_samples/Run3V02/ntuples_xsec_sf_mu_{v}/2022/EGamma_Run2022C-PromptReco-v1/ee/*.root'.format(v=v)]
print(path_mc)

event_selection = "(q_1*q_2 < 0)"
event_selection += " && (trg_single_mu24_1 || trg_single_mu24_2)"
#event_selection += " && (trg_single_ele27_1 || trg_single_ele27_2)"




rdf_mc = ROOT.RDataFrame("ntuple", path_mc)
rdf_dt = ROOT.RDataFrame("ntuple", path_dt)

rdf_mc = rdf_mc.Filter(event_selection)
rdf_dt = rdf_dt.Filter(event_selection)    

#hist_mc = rdf_mc.Histo1D(('mc', "m_vis_corr", 50, 80, 110), 'm_vis_corr_y15')
hist_mc = rdf_mc.Histo1D(('mc', "m_vis_corr", 200, 50, 130), 'm_vis_corr')
#hist_dt = rdf_dt.Histo1D(('data', "m_vis_corr", 50, 80, 110), 'm_vis_corr_y15')
hist_dt = rdf_dt.Histo1D(('data', "m_vis_corr", 200, 50, 130), 'm_vis_corr')

file0 = ROOT.TFile('test.root', 'RECREATE')
hist_mc.Write()
hist_dt.Write()
file0.Close()

file0 = ROOT.TFile('test.root')
hist_mc, hist_dt = file0.Get('mc'), file0.Get('data')

hist_mc.Scale(1./hist_mc.Integral())
evts=hist_dt.Integral()
hist_dt.Scale(1./evts)

plots={'mc': hist_mc, 'dt': hist_dt}

#hist_dt.Divide(hist_mc)
plots['mc'].GetYaxis().SetRangeUser(0,.11)

chi2 = utils.plot_ratio(plots=plots, rcolors={'mc': ROOT.kBlue, 'dt': ROOT.kBlack}, title="total m_vis_corr", outfile="test_ratio_{}.pdf".format(v), evts=evts)


"""
filter_template = "({bin}_{n}{corr} > {bin_l} && {bin}_{n}{corr} < {binr})"
bins = {'eta': [-2.4, -1.2, 0., 1.2, 2.4], 'pt': [25, 30, 35, 39, 42, 45, 49, 60, 1500]}
bin1, bin2 = list(bins.keys())[0], list(bins.keys())[1]
nbins1, nbins2 = len(bins[bin1])-1, len(bins[bin2])-1
c='_corr'


rdf_mc = rdf_mc.Define("diff_pt1", "pt_1_corr - pt_1").Define("diff_pt2", "pt_2_corr - pt_2").Define("diff_m", "m_vis_corr - m_vis")
#rdf_dt = rdf_dt.Define("diff_pt1", "pt_1_corr - pt_1").Define("diff_pt2", "pt_2_corr - pt_2").Define("diff_m", "m_vis_corr - m_vis")

for i in range(nbins1):
    for j in range(nbins2):

        filter1a = filter_template.format(bin=bin1, n=1, bin_l=bins[bin1][i], binr = bins[bin1][i+1], corr='')
        filter1b = filter_template.format(bin=bin2, n=1, bin_l=bins[bin2][j], binr = bins[bin2][j+1], corr=c)
        filter2a = filter_template.format(bin=bin1, n=2, bin_l=bins[bin1][i], binr = bins[bin1][i+1], corr='')
        filter2b = filter_template.format(bin=bin2, n=2, bin_l=bins[bin2][j], binr = bins[bin2][j+1], corr=c)

        rdf_help = rdf_mc.Filter("( {} && {} )".format(filter1a, filter1b))
        mean_p1 = rdf_help.Mean("pt_1_corr").GetValue()
        std_p1  = rdf_help.StdDev("diff_pt1").GetValue()

        rdf_help = rdf_mc.Filter("( {} && {} )".format(filter2a, filter2b))
        mean_p2  = rdf_help.Mean("pt_2_corr").GetValue()
        std_p2   = rdf_help.StdDev("diff_pt2").GetValue()

        rdf_help = rdf_mc.Filter("( {} && {} ) || ( {} && {} )".format(filter1a, filter1b, filter2a, filter2b))
        mean_m = rdf_help.Mean("diff_m").GetValue()
        std_m  = rdf_help.StdDev("diff_m").GetValue()

        smear_true = np.sqrt(mz_res_dt[i][j]**2 - mz_res_mc[i][j]**2)

        print("########################################################")
        print("eta bin: ", i, "; pt bin: ", j)
        print("needed mass smearing: ", smear_true)
        print("Lepton1 - calculated mean lepton smearing: ", mean_p1 * 2 * smear_true/91.1876, " vs true smearing ", std_p1)
        print("Lepton2 - calculated mean lepton smearing: ", mean_p2 * 2 * smear_true/91.1876, " vs true smearing ", std_p2)

        print("projected mass smearing: ", 91.1876/2*(std_p1/mean_p1 + std_p2/mean_p2))
        print("real mass smearing: ", std_m)
        print("########################################################")

#######

"""
"""
c = ROOT.TCanvas("c", "title", 600, 600)
#pt_mc = rdf_mc.Filter("({} && {})".format(filter2a, filter2b)).Histo1D(("p2_mc", "pt_mc", 50, bins['pt'][b2]-2, bins['pt'][b2+1]+2), "m_vis_corr")
#pt_dt = rdf_dt.Filter("({} && {})".format(filter2a, filter2b)).Histo1D(("p2_dt", "pt_dt", 50, bins['pt'][b2]-2, bins['pt'][b2+1]+2), "m_vis_corr")
pt_mc = rdf_mc.Filter("({} && {}) || ({} && {})".format(filter1a, filter1b, filter2a, filter2b)).Histo1D(("p2_mc", "pt_mc", 200, 50, 130), "m_vis_corr")
pt_dt = rdf_dt.Filter("({} && {}) || ({} && {})".format(filter1a, filter1b, filter2a, filter2b)).Histo1D(("p2_dt", "pt_dt", 200, 50, 130), "m_vis_corr")


pt_mc.Scale(1./pt_mc.Integral())
pt_dt.Scale(1./pt_dt.Integral())
pt_dt.SetMarkerStyle(20)
pt_mc.Draw("HIST")
pt_mc.GetYaxis().SetRangeUser(10e-5, .1)
pt_dt.Draw("same")



c.SetLogy(1)
c.SaveAs("test.pdf")
"""