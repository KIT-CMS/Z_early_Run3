#from winreg import HKEY_DYN_DATA
import ROOT
import os,sys
import numpy as np
from CMSPLOTS.myFunction import DrawHistos

doMuon = True
doElectron = False
doWpT = False

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(ROOT.kFALSE)
ROOT.TH2.AddDirectory(ROOT.kFALSE)
ROOT.TH3.AddDirectory(ROOT.kFALSE)

etaLabels = {
    "": "",
    "barrel": "Barrel", 
    "endcap": "Endcap"
}
'''
channelLabels = {
    "muplus":  "W^{+}#rightarrow #mu^{+}#nu",
    "muminus": "W^{-}#rightarrow #mu^{-}#nu",
    "eplus":   "W^{+}#rightarrow e^{+}#nu",
    "eminus":  "W^{-}#rightarrow e^{-}#nu",
}
'''
if doWpT:
    wptbins = ["WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8", "WpT_bin9"]
else:
    wptbins = ["WpT_bin0"]

    #TODO: MORE THAN MT VARIABLE


def ExpltOneBin(variable, isocenters, bincontents, binerrors, isoSR, mTmin, mTmax, suffix="", extraText="", bincontents_scaled = None, binerrors_scaled = None):
    """
    extrapolate the QCD shape from a set of control regions (isocenters) to the signal region (isoSR),
    using linear extrapolation and the 2nd order polynomial function
    """
    graph = ROOT.TGraphErrors(len(bincontents), np.array(isocenters), np.array(bincontents), np.zeros(len(bincontents)), np.array(binerrors))
    graph.SetMarkerStyle(20)
    f0 = ROOT.TF1("pol0_"+suffix, "[0]", -0.1, 1.0)
    f1 = ROOT.TF1("pol1_"+suffix, "[0]*(x-{})+[1]".format(str(isoSR)), -0.1, 1.0)
    f2 = ROOT.TF1("pol2_"+suffix, "[0]*(x-{isoSR})*(x-{isoSR})+[1]*(x-{isoSR})+[2]".format(isoSR=str(isoSR)), -0.1, 1.0)
    # fit range
    fitmin = 0.20  # 0.25
    fitmax = 1.00  # 0.60
    graph.Fit(f0, "R", "", fitmin, fitmax)
    graph.Fit(f1, "R", "", fitmin, fitmax)
    # graph.Fit(f2, "R", "", fitmin, fitmax)
    #print("val at Signal region", f1.Eval(isoSR))

    val_pol0_par1 = f0.GetParameter(0)
    err_pol0_par1 = f0.GetParError(0)
    val_pol1_par1 = f1.GetParameter(1)
    err_pol1_par1 = f1.GetParError(1)
    val_pol1_par0 = f1.GetParameter(0)
    err_pol1_par0 = f1.GetParError(0)
    # val_pol2_par2 = f2.GetParameter(2)
    # err_pol2_par2 = f2.GetParError(2)
    #print("val ", val, " error ", err)

    f0.SetLineStyle(2)
    f0.SetLineColor(25)
    f1.SetLineStyle(2)
    f1.SetLineColor(46)
    # f2.SetLineStyle(2)
    # f2.SetLineColor(9)

    graph1 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val_pol0_par1]), np.zeros(1), np.array([abs(err_pol0_par1)]))
    graph1.SetMarkerColor(25)
    graph1.SetMarkerStyle(47)
    graph1.SetMarkerSize(2)

    graph2 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val_pol1_par1]), np.zeros(1), np.array([abs(err_pol1_par1)]))
    graph2.SetMarkerColor(46)
    graph2.SetMarkerStyle(47)
    graph2.SetMarkerSize(2)

    # graph3 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val_pol2_par2]), np.zeros(1), np.array([abs(err_pol2_par2)]))
    # graph3.SetMarkerColor(9)
    # graph3.SetMarkerStyle(45)
    # graph3.SetMarkerSize(2)

    bin_str = None
    if "ptOverPfmt" in variable:
        bin_str = "{} < pT/mT < {}".format(mTmin, mTmax)
    elif "mt" in variable:
        bin_str = "{} < mT < {}".format(mTmin, mTmax)
    elif "met" in variable:
        bin_str = "{} < MET < {}".format(mTmin, mTmax)
    h_todraws = [graph, f0, f1, graph1, graph2]
    labels = [bin_str, "Pol0 Fit", "Pol1 Fit", "Pol0 Extrapolation", "Pol1 Extrapolation"]
    drawoptions = ["P same", "L", "L", "P same", "P same"]
    legendoptions=["EP", "L", "L", "EP", "EP"]
    # h_todraws = [graph, f0, f1, f2, graph1, graph2, graph3]
    # labels = [bin_str, "Pol0 Fit", "Pol1 Fit", "Pol2 Fit", "Pol0 Extrapolation", "Pol1 Extrapolation", "Pol2 Extrapolation"]
    # drawoptions = ["P same", "L", "L", "L", "P same", "P same", "P same"]
    # legendoptions=["EP", "L", "L", "L", "EP", "EP", "EP"]

    val_scaled_pol1_par1 = None
    err_scaled_pol1_par1 = None
    if bincontents_scaled:
        # repeat the fitting precedure on the scaled templates
        graph_scaled = ROOT.TGraphErrors(len(bincontents), np.array(isocenters), np.array(bincontents_scaled), np.zeros(len(bincontents)), np.array(binerrors_scaled))
        f3 = ROOT.TF1("pol1_scaled"+suffix, "[0]*(x-{})+[1]".format(str(isoSR)), -0.1, 0.60)
        graph_scaled.Fit(f3, "R", "", fitmin, fitmax)
        val_scaled_pol1_par1 = f3.GetParameter(1)
        err_scaled_pol1_par1 = f3.GetParError(1)

        f3.SetLineStyle(2)
        f3.SetLineColor(30)

        graph_scaled.SetMarkerColor(12)
        graph_scaled.SetMarkerStyle(31)
        graph4 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val_scaled_pol1_par1]), np.zeros(1), np.array([abs(err_scaled_pol1_par1)]))
        graph4.SetMarkerColor(30)
        graph4.SetMarkerStyle(41)
        graph4.SetMarkerSize(2)

        h_todraws.append(graph_scaled)
        h_todraws.append(f3)
        h_todraws.append(graph4)

        labels.append("Templates with Scaled MC")
        labels.append("Pol1 Fit with Scaled MC")
        labels.append("Pol1 Extrapolation with Scaled MC")
        
        drawoptions.extend(["P same", "L", "P same"])
        legendoptions.extend(["PE", "L", "PE"])

    bincontents_forAxis = bincontents
    bincontents_forAxis.append(val_pol0_par1)
    bincontents_forAxis.append(val_pol1_par1)
    # bincontents_forAxis.append(val_pol2_par2)    

    DrawHistos(
        h_todraws, labels, 0, 1.6, "Lepton Relative Isolation",
        0.2*min(bincontents_forAxis),
        1.5*max(bincontents_forAxis),
        "Bin Content", "QCDBinContentNorm_"+suffix,
        dology=False,
        drawoptions=drawoptions,
        legendoptions=legendoptions,
        nMaxDigits=3,
        legendPos=[0.65, 0.18, 0.88, 0.58],
        lheader=extraText
    )

    #return (val_pol1_par1, err_pol1_par1), (val_pol1_par0, err_pol1_par0), (val_pol2_par2, err_pol2_par2), (val_scaled_pol1_par1, err_scaled_pol1_par1)
    return (val_pol1_par1, err_pol1_par1), (val_pol1_par0, err_pol1_par0), (val_pol0_par1, err_pol0_par1)


def ExtrapolateQCD(fname, oname, channel, variable, wptbin, etabins, isobins, mc_scale = 1.):
    """
    run the QCd extrapolation in all mT bins,
    save the statistical and systematic variations for HComb, and
    make some shape comparison plots
    """

    assert mc_scale in [1., 1.3, 0.7]  # nominal up / down
    mc_scale_str = ""
    if mc_scale > 1.:
        mc_scale_str = "_mcScaleUp"
    elif mc_scale < 1.:
        mc_scale_str = "_mcScaleDown"

    if not os.path.exists("root/QCD"):
        os.makedirs("root/QCD")
    ofile = ROOT.TFile.Open("root/QCD/"+oname, "recreate")

    #ADDED BY JULIUS
    #examples of bins: [w_iso6], [lepEta_bin0], [WpT_bin0]

    isomin = isobins[0]
    isomax = isobins[-1]
    #TODO: Make sure these are right
    isocuts = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]
    isocenters = [(isocuts[i]+isocuts[i+1])/2 for i in range(len(isocuts)-1)]
    
    isoSR = 0.032  # 0.025 # average isolation value in the signal region

    h_isoSR = None
    h_isoSR_noMCscale = None
    for etabin in etabins:
        histos_norm = dict()
        histos_scaled_norm = dict()
        fqcd = ROOT.TFile.Open(fname)
        fqcd.Print()

        for iso in ["isoSR"]+isobins:
            if doMuon:
                hname = "data#mmet_corr-355862-357482-data#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                DY_hname = "DY#mmet_corr-355862-357482-DY#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                W_hname = "W#mmet_corr-355862-357482-W#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                TT_hname = "TT#mmet_corr-355862-357482-TT#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                ST_hname = "ST#mmet_corr-355862-357482-ST#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                VV_hname = "VV#mmet_corr-355862-357482-VV#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                DYtau_hname = "DYtau#mmet_corr-355862-357482-DYtau#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                Wtau_hname = "Wtau#mmet_corr-355862-357482-Wtau#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                print(hname)
                h = fqcd.Get(hname)
                h_DY = fqcd.Get(DY_hname)
                h_W = fqcd.Get(W_hname)
                h_TT = fqcd.Get(TT_hname)
                h_ST = fqcd.Get(ST_hname)
                h_VV = fqcd.Get(VV_hname)
                h_DYtau = fqcd.Get(DYtau_hname)
                h_Wtau = fqcd.Get(Wtau_hname)

                #Combine all MC histograms
                h_mc_total_noMCscale = h_DY.Clone("h_mc_total_noMCscale")
                h_mc_total_noMCscale.Add(h_W)
                h_mc_total_noMCscale.Add(h_TT)
                h_mc_total_noMCscale.Add(h_ST)
                h_mc_total_noMCscale.Add(h_VV)
                h_mc_total_noMCscale.Add(h_DYtau)
                h_mc_total_noMCscale.Add(h_Wtau)

                # scale signal MC
                h_W.Scale(mc_scale)

                h_mc_total = h_DY.Clone("h_mc_total")
                h_mc_total.Add(h_W)
                h_mc_total.Add(h_TT)
                h_mc_total.Add(h_ST)
                h_mc_total.Add(h_VV)
                h_mc_total.Add(h_DYtau)
                h_mc_total.Add(h_Wtau)

                # remove negative entries from total MC histogram
                for ibin in range(0, h_mc_total.GetNbinsX()+2):
                    h_mc_total_noMCscale.SetBinContent(ibin, max(h_mc_total_noMCscale.GetBinContent(ibin), 0.))
                    h_mc_total.SetBinContent(ibin, max(h_mc_total.GetBinContent(ibin), 0.))

                if iso == "isoSR":
                    h_isoSR_noMCscale = h.Clone("h_isoSR_noMCscale")
                    h_isoSR_noMCscale.Add(h_mc_total_noMCscale, -1.)

                # subtract MC contribution from data histogram
                h.Add(h_mc_total, -1)
                
                # remove negative entries after the subtraction
                for ibin in range(0, h.GetNbinsX()+2):
                    h.SetBinContent(ibin, max(h.GetBinContent(ibin), 0.))
                    if iso == "isoSR":
                        h_isoSR_noMCscale.SetBinContent(ibin, max(h_isoSR_noMCscale.GetBinContent(ibin), 0.))
                h.Print()
            elif doElectron:
                hname = "data#emet#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                DY_hname = "DY#emet-DY#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                W_hname = "W#emet-W#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                TT_hname = "TT#emet-TT#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                ST_hname = "ST#emet-ST#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                VV_hname = "VV#emet-VV#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                DYtau_hname = "DYtau#emet-DYtau#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                Wtau_hname = "Wtau#emet-Wtau#Nominal#{}_{}{}_{}".format(variable, channel, etabin, iso)
                print(hname)
                h = fqcd.Get(hname)
                h_DY = fqcd.Get(DY_hname)
                h_W = fqcd.Get(W_hname)
                h_TT = fqcd.Get(TT_hname)
                h_ST = fqcd.Get(ST_hname)
                h_VV = fqcd.Get(VV_hname)
                h_DYtau = fqcd.Get(DYtau_hname)
                h_Wtau = fqcd.Get(Wtau_hname)

                #Combine all MC histograms
                h_mc_total_noMCscale = h_DY.Clone("h_mc_total_noMCscale")
                h_mc_total_noMCscale.Add(h_W)
                h_mc_total_noMCscale.Add(h_TT)
                h_mc_total_noMCscale.Add(h_ST)
                h_mc_total_noMCscale.Add(h_VV)
                h_mc_total_noMCscale.Add(h_DYtau)
                h_mc_total_noMCscale.Add(h_Wtau)

                # scale signal MC
                h_W.Scale(mc_scale)

                h_mc_total = h_DY.Clone("h_mc_total")
                h_mc_total.Add(h_W)
                h_mc_total.Add(h_TT)
                h_mc_total.Add(h_ST)
                h_mc_total.Add(h_VV)
                h_mc_total.Add(h_DYtau)
                h_mc_total.Add(h_Wtau)

                # remove negative entries from total MC histogram
                for ibin in range(0, h_mc_total.GetNbinsX()+2):
                    h_mc_total_noMCscale.SetBinContent(ibin, max(h_mc_total_noMCscale.GetBinContent(ibin), 0.))
                    h_mc_total.SetBinContent(ibin, max(h_mc_total.GetBinContent(ibin), 0.))

                if iso == "isoSR":
                    h_isoSR_noMCscale = h.Clone("h_isoSR_noMCscale")
                    h_isoSR_noMCscale.Add(h_mc_total_noMCscale, -1.)

                # subtract MC contribution from data histogram
                h.Add(h_mc_total, -1)
                
                # remove negative entries after the subtraction
                for ibin in range(0, h.GetNbinsX()+2):
                    h.SetBinContent(ibin, max(h.GetBinContent(ibin), 0.))
                    if iso == "isoSR":
                        h_isoSR_noMCscale.SetBinContent(ibin, max(h_isoSR_noMCscale.GetBinContent(ibin), 0.))
                h.Print()

            # set the overflow and underflow to zero
            h.SetBinContent(0, 0)
            h.SetBinContent(h.GetNbinsX()+1, 0)
            h_norm =  h.Clone("h_norm_"+channel+"_"+iso)
            if iso == "isoSR":
                h_isoSR = h_norm
            else:
                histos_norm[iso] = h_norm

        href = h_isoSR_noMCscale  # h_data_isoSR  # h_isoSR  # histos_norm["isoSR"]
        counts = href.Integral()
        for iso in isobins:
            #Normalize histrograms
            # histos_norm[iso].Scale(1. / histos_norm[iso].Integral())
            histos_norm[iso].Scale(counts / histos_norm[iso].Integral())

        # some histograms to save the trend of function parameter variation as a function of mT
        # mostly for plotting purpose
        h_pol0_par1 = href.Clone("h_pol0_par1_{}_{}_{}".format(channel, etabin, wptbin)) # interception
        h_pol1_par1 = href.Clone("h_pol1_par1_{}_{}_{}".format(channel, etabin, wptbin)) # interception
        h_pol1_par0 = href.Clone("h_pol1_par0_{}_{}_{}".format(channel, etabin, wptbin)) # slope
        h_pol2_par2 = href.Clone("h_pol2_par2_{}_{}_{}".format(channel, etabin, wptbin)) # interception
        #h_scaled_pol1_par1 = href.Clone("h_scaled_pol1_par1_{}_{}_{}".format(channel, etabin, wptbin)) # interception for the scaled templates

        # save the extrapolated shape for HComb
        hnew = href.Clone("h_QCD_Extrapolated_"+variable+"_"+channel)
        hnew_pol1 = href.Clone("h_QCD_Extrapolated_"+variable+"_"+channel+"_Pol1shapeUp")
        hnew_pol2 = href.Clone("h_QCD_Extrapolated_"+variable+"_"+channel+"_Pol2shapeUp")
        hnew_scaled = href.Clone("h_QCD_Extrapolated_"+variable+"_"+channel+"_ScaledMCshapeUp")


        vals_pol0_par1 = []
        vals_pol1_par1 = []
        #
        # run the linear extrapolation bin-by-bin
        #
        for ibin in xrange(1, histos_norm[isomin].GetNbinsX()+1):
            bincontents = []
            binerrors = []
            for iso, hist in histos_norm.iteritems():
                bincontents.append( hist.GetBinContent(ibin) )
                binerrors.append( hist.GetBinError(ibin) )

            mTmin = histos_norm[isomin].GetBinLowEdge(ibin)
            mTmax = histos_norm[isomin].GetBinLowEdge(ibin)+histos_norm[isomin].GetBinWidth(ibin)

            suffix = variable+"_"+channel+"_bin_"+str(ibin)+mc_scale_str  # +"_"+ etabin+"_"+wptbin
            if doMuon and channel == "pos":
                extraText = "W^{+}#rightarrow #mu^{+}#nu" +" "+ etaLabels[etabin]
            elif doMuon and channel == "neg":
                extraText = "W^{-}#rightarrow #mu^{-}#nu" +" "+ etaLabels[etabin]
            elif not doMuon and channel == "pos":
                extraText = "W^{+}#rightarrow e^{+}#nu" +" "+ etaLabels[etabin]
            elif not doMuon and channel == "neg":
                extraText = "W^{-}#rightarrow e^{-}#nu" +" "+ etaLabels[etabin]
 

            bincontents_scaled = None
            binerrors_scaled = None
            results_pol1_par1, results_pol1_par0, results_pol0_par1 = ExpltOneBin(variable, isocenters, bincontents, binerrors, isoSR, mTmin, mTmax, suffix=suffix, extraText=extraText, bincontents_scaled = bincontents_scaled, binerrors_scaled = binerrors_scaled)

            hnew.SetBinContent(ibin, max(results_pol0_par1[0], 0))
            hnew.SetBinError(ibin, 0.)
            # hnew.SetBinContent(ibin, max(results_pol1_par1[0], 0))
            # hnew.SetBinError(ibin, 0.)

            hnew_pol1.SetBinContent(ibin, max(results_pol1_par1[0], 0))
            hnew_pol1.SetBinError(ibin, 0.)

            h_pol0_par1.SetBinContent(ibin, results_pol0_par1[0])
            h_pol0_par1.SetBinError(ibin,   results_pol0_par1[1])
            h_pol1_par1.SetBinContent(ibin, results_pol1_par1[0])
            h_pol1_par1.SetBinError(ibin,   results_pol1_par1[1])
            h_pol1_par0.SetBinContent(ibin, results_pol1_par0[0])
            h_pol1_par0.SetBinError(ibin,   results_pol1_par0[1])
            # h_pol2_par2.SetBinContent(ibin, results_pol2_par2[0])
            # h_pol2_par2.SetBinError(ibin,   results_pol2_par2[1])

            vals_pol0_par1.append(results_pol0_par1)
            vals_pol1_par1.append(results_pol1_par1)

        # set the bin-by-bin shape variation (stat.) for HComb
        hnew_ups = []
        hnew_downs = []
        for ibin in xrange(1, histos_norm[isomin].GetNbinsX()+1):
            val = max(vals_pol0_par1[ibin-1][0], 0.)
            stat = vals_pol0_par1[ibin-1][1]
            pol1 = hnew.GetBinContent(ibin) - hnew_pol1.GetBinContent(ibin)
            err = np.sqrt(stat**2 + pol1**2)
            hnew_up   = hnew.Clone("h_QCD_Extrapolated_"+variable+"_"+channel+"_bin{}shapeUp".format(str(ibin)))
            hnew_down = hnew.Clone("h_QCD_Extrapolated_"+variable+"_"+channel+"_bin{}shapeDown".format(str(ibin)))
            hnew_up.SetBinContent(ibin, val+err)
            hnew_up.SetBinError(ibin, 0.)
            hnew_down.SetBinContent(ibin, max(val-err, 0.))
            hnew_down.SetBinError(ibin, 0.)

            hnew_ups.append(hnew_up)
            hnew_downs.append(hnew_down)

        #STOPPED ADDED BY JULIUS

        # pol1 as another systematic
        hnew_pol1.Scale(hnew.Integral() / hnew_pol1.Integral())
        hnew_pol1Dn = hnew_pol1.Clone("h_QCD_Extrapolated_"+variable+"_"+channel+"_Pol1shapeDown")
        for ibin in xrange(1, hnew.GetNbinsX()+1):
            hnew_pol1Dn.SetBinContent(ibin, max(0., 2*hnew.GetBinContent(ibin) - hnew_pol1.GetBinContent(ibin)))
        hnew_ups.append(hnew_pol1)
        hnew_downs.append(hnew_pol1Dn)

        # pol2 as another systematic
        # hnew_pol2.Scale(hnew.Integral() / hnew_pol2.Integral())
        # hnew_pol2Dn = hnew_pol2.Clone("h_QCD_Extrapolated_"+variable+"_"+channel+"_Pol2shapeDown")
        # for ibin in xrange(1, hnew.GetNbinsX()+1):
        #     hnew_pol2Dn.SetBinContent(ibin, 2*hnew.GetBinContent(ibin) - hnew_pol2.GetBinContent(ibin))
        # hnew_ups.append(hnew_pol2)
        # hnew_downs.append(hnew_pol2Dn)

        # scaled MC as another systematic
        # if fname_scaled:
        #     hnew_scaled.Scale(hnew.Integral() / hnew_scaled.Integral())
        #     hnew_scaledDn = hnew_scaled.Clone("h_QCD_Extrapolated_"+variable+"_"+channel+"_ScaledMCshapeDown")
        #     for ibin in xrange(1, hnew.GetNbinsX()+1):
        #         hnew_scaledDn.SetBinContent(ibin, max(0., 2*hnew.GetBinContent(ibin) - hnew_scaled.GetBinContent(ibin)))
        #     hnew_ups.append(hnew_scaled)
        #     hnew_downs.append(hnew_scaledDn)

        h_pol0_par1.SetLineColor(25)
        h_pol0_par1.SetMarkerColor(25)
        h_pol1_par1.SetLineColor(46)
        h_pol1_par1.SetMarkerColor(46)
        h_pol1_par1.Scale(h_pol0_par1.Integral() / h_pol1_par1.Integral())
        h_pol2_par2.Scale(h_pol0_par1.Integral() / h_pol2_par2.Integral())
        h_pol2_par2.SetLineColor(9)
        h_pol2_par2.SetMarkerColor(9)
        h_todraws = [h_pol0_par1, h_pol1_par1]  # , h_pol2_par2]
        labels = ["Pol0 Extrapolation", "Pol1 Extrapolation"]  #, "Pol2 Extrapolation"]
        # if fname_scaled:
        #     h_scaled_pol1_par1.SetLineColor(30)
        #     h_scaled_pol1_par1.SetMarkerColor(30)
        #     h_todraws.append(h_scaled_pol1_par1)
        #     labels.append("Scaled MC")

        if "ptOverPfmt" in variable:
            DrawHistos( h_todraws, labels, 0, 3, "#frac{p_{T}}{m_{T}}", 0., 1.25*h_isoSR_noMCscale.GetMaximum(), "A.U.", "QCDShapeCompare_"+variable+"_"+channel+mc_scale_str, dology=False, nMaxDigits=3, legendPos=[0.60, 0.72, 0.88, 0.88], lheader=extraText, noCMS = True)
        elif "mt" in variable:
            #DrawHistos( h_todraws, labels, 0, 120, "m_{T} [GeV]", 0., 1.25*h_isoSR_noMCscale.GetMaximum(), "A.U.", "QCDShapeCompare_"+variable+"_" +channel+mc_scale_str, dology=False, nMaxDigits=3, legendPos=[0.60, 0.72, 0.88, 0.88], lheader=extraText, noCMS = True)
            DrawHistos( h_todraws, labels, 0, 120, "m_{T} [GeV]", 0., 1.25*h_isoSR_noMCscale.GetMaximum(), "A.U.", "QCDShapeCompare_"+variable+"_" +channel+mc_scale_str, dology=False, nMaxDigits=3, legendPos=[0.60, 0.72, 0.88, 0.88], lheader=extraText, noCMS = True)
        elif "met" in variable:
            DrawHistos( h_todraws, labels, 0, 120, "MET [GeV]", 0., 1.25*h_isoSR_noMCscale.GetMaximum(), "A.U.", "QCDShapeCompare_"+variable+"_"+channel+mc_scale_str, dology=False, nMaxDigits=3, legendPos=[0.60, 0.72, 0.88, 0.88], lheader=extraText, noCMS = True)

        # 
        # draw the original (normalized) mT distribution in different anti-isolated regions
        #
        h_todraws = []
        labels = []
        suffix = channel
        ROOT.gStyle.SetPalette(77)
        colors = ROOT.TColor.GetPalette()
        for idx, (iso, hmt) in enumerate(histos_norm.items()):
            # if idx >= 9:
            #     continue
            icolor = int(idx*(255./len(isobins)))
            hmt.SetLineColor(colors.At(icolor))
            h_todraws.append(hmt)
            labels.append("{:.2f} < iso < {:.2f}".format(isocuts[idx], isocuts[idx+1]))

        if "ptOverPfmt" in variable:
            DrawHistos( h_todraws, labels, 0, 3, "#frac{p_{T}}{m_{T}}", 0., 0.1, "A.U.", "shapes_QCD_"+variable+"_"+suffix, dology=False, nMaxDigits=3, legendPos=[0.75, 0.40, 0.95, 0.83], lheader=extraText, noCMS = True, donormalize=True, drawashist=True, showratio=True, ratiobase=-1, yrlabel="#scale[0.4]{Ratio to last iso bin}", yrmin=0.51, addOverflow = True, yrmax=1.49)
        elif "mt" in variable:
            #DrawHistos( h_todraws, labels, 0, 120, "m_{T} [GeV]", 0., 0.1, "A.U.", "shapes_QCD_"+variable+"_"+suffix, dology=False, nMaxDigits=3, legendPos=[0.75, 0.40, 0.95, 0.83], lheader=extraText, noCMS = True, donormalize=True, drawashist=True, showratio=True, ratiobase=-1, yrlabel="#scale[0.4]{Ratio to last iso bin}", yrmin=0.51, addOverflow = True, yrmax=1.49)
            DrawHistos( h_todraws, labels, 0, 120, "m_{T} [GeV]", 0., 0.1, "A.U.", "shapes_QCD_"+variable+"_"+suffix, dology=False, nMaxDigits=3, legendPos=[0.75, 0.40, 0.95, 0.83], lheader=extraText, noCMS = True, donormalize=True, drawashist=True, showratio=True, ratiobase=-1, yrlabel="#scale[0.4]{Ratio to last iso bin}", yrmin=0.51, addOverflow = True, yrmax=1.49)
        elif "met" in variable:
            DrawHistos( h_todraws, labels, 0, 120, "MET [GeV]", 0., 0.1, "A.U.", "shapes_QCD_"+variable+"_"+suffix, dology=False, nMaxDigits=3, legendPos=[0.75, 0.40, 0.95, 0.83], lheader=extraText, noCMS = True, donormalize=True, drawashist=True, showratio=True, ratiobase=-1, yrlabel="#scale[0.4]{Ratio to last iso bin}", yrmin=0.51, addOverflow = True, yrmax=1.49)


        #
        # write the variations to the output
        #
        ofile.cd()
        hnew.SetDirectory(ofile)
        hnew.Write(hnew.GetName()+mc_scale_str)
        # for h in h_todraws+[h_pol1_par0]:
        #     h.SetDirectory(ofile)
        #     h.Write(h.GetName()+mc_scale_str)
        for hnew_up in hnew_ups:
            hnew_up.SetDirectory(ofile)
            hnew_up.Write(hnew_up.GetName()+mc_scale_str)
        for hnew_down in hnew_downs:
            hnew_down.SetDirectory(ofile)
            hnew_down.Write(hnew_down.GetName()+mc_scale_str)

        # for iso in isobins:
        #     histos_norm[iso].SetDirectory(ofile)
        #     histos_norm[iso].Write(histos_norm[iso].GetName()+mc_scale_str)

        number_of_qcd_events = h_todraws[0].Integral()
    
    ofile.Close()
    return number_of_qcd_events


if __name__ == "__main__":
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.gROOT.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);")
    ROOT.gROOT.ProcessLine("RooMsgService::instance().setSilentMode(true);")

    isobins = ["iso5", "iso6", "iso7", "iso8", "iso9", "iso10", "iso11", "iso12", "iso13", "iso14", "iso15", "iso16", "iso17", "iso18", "iso19", "iso20"]
    lepEtaBins = [""]  # ["_barrel", "_endcap"]
    met_tag = "pf"
    variables = [met_tag+"mt_corr"] # , met_tag+"met_corr"]  # , "ptOverPfmt_corr"

    number_of_qcd_events_pos = dict()
    number_of_qcd_events_neg = dict()
    number_of_qcd_events_normed = dict()

    if doMuon:
        fname = "/work/jdriesch/earlyrun3/Z_early_Run3/output/earlyRun3_2022_QCD.root"
        for variable in variables:
            for wptbin in ["WpT_bin0"]:
                oname = "qcdshape_extrapolated_muplus_"+variable
                events = ExtrapolateQCD(fname, oname+".root", "pos", variable, wptbin, lepEtaBins, isobins, 1.)
                number_of_qcd_events_pos[variable] = events
                events = ExtrapolateQCD(fname, oname+"_mcScaleUp.root", "pos", variable, wptbin, lepEtaBins, isobins, 1.3)
                events = ExtrapolateQCD(fname, oname+"_mcScaleDown.root", "pos", variable, wptbin, lepEtaBins, isobins, 0.7)

                oname = "qcdshape_extrapolated_muminus_"+variable
                events = ExtrapolateQCD(fname, oname+".root", "neg", variable, wptbin, lepEtaBins, isobins, 1.)
                number_of_qcd_events_neg[variable] = events
                events = ExtrapolateQCD(fname, oname+"_mcScaleUp.root", "neg", variable, wptbin, lepEtaBins, isobins, 1.3)
                events = ExtrapolateQCD(fname, oname+"_mcScaleDown.root", "neg", variable, wptbin, lepEtaBins, isobins, 0.7)


        
        # for variable in variables:
        #     number_of_qcd_events_normed["{}_pos".format(variable)] = number_of_qcd_events_pos[variable] / number_of_qcd_events_pos[met_tag+"met_corr"]
        #     number_of_qcd_events_normed["{}_neg".format(variable)] = number_of_qcd_events_neg[variable] / number_of_qcd_events_pos[met_tag+"met_corr"]
        # print("Positive: ", number_of_qcd_events_pos)
        # print("Negative: ", number_of_qcd_events_neg)
        # print(number_of_qcd_events_normed)


    # if doElectron:
        # fname = "XXX"
        # for variable in variables:
            # for wptbin in ["WpT_bin0"]:
                # oname = "qcdshape_extrapolated_eplus_"+variable
                # events = ExtrapolateQCD(fname, oname+"_"+wptbin+"_eplus.root", "pos", variable, wptbin, lepEtaBins, isobins)
                # number_of_qcd_events_pos[variable] = events

                # oname = "qcdshape_extrapolated_eminus_"+variable
                # events = ExtrapolateQCD(fname, oname+"_"+wptbin+"_eminus.root", "neg", variable, wptbin, lepEtaBins, isobins)
                # number_of_qcd_events_neg[variable] = events

        # for variable in variables:
            # number_of_qcd_events_normed["{}_pos".format(variable)] = number_of_qcd_events_pos[variable] / number_of_qcd_events_pos[met_tag+"met_corr"]
            # number_of_qcd_events_normed["{}_neg".format(variable)] = number_of_qcd_events_neg[variable] / number_of_qcd_events_pos[met_tag+"met_corr"]
        # print("Positive: ", number_of_qcd_events_pos)
        # print("Negative: ", number_of_qcd_events_neg)
        # print(number_of_qcd_events_normed)
