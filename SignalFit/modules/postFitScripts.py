"""
script to make postfit comparisons
"""
import ROOT
import os,sys,math
from collections import OrderedDict
import json
import re
import numpy as np

from CMSPLOTS.myFunction import DrawHistos, DrawConfig

ROOT.gROOT.SetBatch(True)

def MakePostPlot(ifilename: str, channel: str, prepost: str, bins: np.array, suffix: str, showpull: bool = False, x_label: str = "", is5TeV: bool = False, startbin: int = 1):
    """
    compare the unrolled postfit of data and templates
    """

    print("")
    print("#"*50)
    print("channel:", channel)
    print("prepost:", prepost)
    print("suffix: ", suffix)
    print("ifile:  ", ifilename)

    ifile = ROOT.TFile(ifilename)

    horgdata = ifile.Get("obs")

    # get the list of histograms saved in the file
    hkeys = ifile.GetListOfKeys()
    hkeys = [hkey.GetName() for hkey in hkeys]
    hnames_sig = []
    hnames_sig_z = []
    hnames_qcd = []
    for hkey in hkeys:
        hkey_str = str(hkey)
        if not hkey_str.endswith(prepost):
            continue

        # w signal
        if bool(re.match(r"expproc_\w*plus_sig_\w*fit$", hkey)) or bool(re.match(r"expproc_\w*minus_sig_\w*fit$", hkey)):
            hnames_sig.append( hkey )
        # z signal
        elif bool(re.match(r"expproc_\w*sig_\w*fit$", hkey)):
            hnames_sig_z.append( hkey )
        # qcd
        elif bool(re.match(r"expproc_qcd_\w*fit$", hkey)):
            hnames_qcd.append( hkey )
    assert len(hnames_sig)>=1, "There should be at least one sig histogram in file: {}".format(ifilename)
    print(f"W signals: {hnames_sig}")
    print(f"Z signals: {hnames_sig_z}")
    print(f"QCD:       {hnames_qcd}")

    # ewk bkg includes W->tau+nu, z->ll, and diboson process (for w's)
    # ewk processes for z's
    hnames_ewks = [f"expproc_ewk_{prepost}"]
    hnames_ttbar = [f"expproc_tt_{prepost}"]

    ## read the postfit plots from input file
    hexpsig = None
    hexpsig_z = None
    hexpewk = None
    hexpqcd = None
    hexpttbar = None

    for hkey in hkeys:
        if hkey in hnames_sig:
            if hexpsig is None:
                hexpsig = ifile.Get(hkey)
            else:
                hexpsig.Add( ifile.Get(hkey) )

        if hkey in hnames_sig_z:
            if hexpsig_z is None:
                hexpsig_z = ifile.Get(hkey)
            else:
                hexpsig_z.Add( ifile.Get(hkey) )

        if hkey in hnames_ewks:
            if hexpewk is None:
                hexpewk = ifile.Get(hkey)
            else:
                hexpewk.Add( ifile.Get(hkey) )

        if hkey in hnames_ttbar:
            if hexpttbar is None:
                hexpttbar = ifile.Get(hkey)
            else:
                hexpttbar.Add( ifile.Get(hkey) )

        if hkey in hnames_qcd:
            if hexpqcd is None:
                hexpqcd = ifile.Get(hkey)
            else:
                hexpqcd.Add( ifile.Get(hkey) )

    # the combined prediction of all processes,
    # which should have included the correct total postfit uncertainties
    hexpfull = ifile.Get(f"expfull_{prepost}")

    # the histograms saved in the root file does not follow the original bining
    # recover the original binning
    nbins = len(bins) - 1
    binnings = (nbins, bins)
    bin_width = bins[1] - bins[0]
    for ibin in range(nbins):
        assert bins[ibin+1] - bins[ibin] == bin_width

    #binnings = (newbins.shape[0]-1, newbins)

    hdata   = ROOT.TH1D("hdata_{}_{}".format( channel, suffix),  "hdata_{}_{}".format( channel, suffix),  *binnings)
    hsig    = ROOT.TH1D("hsig_{}_{}".format(  channel, suffix),  "hsig_{}_{}".format(  channel, suffix),  *binnings)
    hsig_z  = ROOT.TH1D("hsig_z_{}_{}".format(channel, suffix),  "hsig_z_{}_{}".format(channel, suffix),  *binnings)
    hewk    = ROOT.TH1D("hewk_{}_{}".format(  channel, suffix),  "hewk_{}_{}".format(  channel, suffix),  *binnings)
    httbar  = ROOT.TH1D("httbar_{}_{}".format(channel, suffix),  "httbar_{}_{}".format(channel, suffix),  *binnings)
    hqcd    = ROOT.TH1D("hqcd_{}_{}".format(  channel, suffix),  "hqcd_{}_{}".format(  channel, suffix),  *binnings)
    hratio  = ROOT.TH1D("hrato_{}_{}".format( channel, suffix),  "hratio_{}_{}".format(channel, suffix),  *binnings)
    hpull   = ROOT.TH1D("hpull_{}_{}".format( channel, suffix),  "hpull_{}_{}".format( channel, suffix),  *binnings)
    for ibin in range(1, nbins + 1):
        hdata.SetBinContent(ibin, horgdata.GetBinContent(ibin + startbin-1))
        hdata.SetBinError(ibin,   horgdata.GetBinError(ibin + startbin-1 ))
        if hexpsig:
            hsig.SetBinContent(ibin,    hexpsig.GetBinContent(ibin + startbin-1))
        if hexpsig_z:
            hsig_z.SetBinContent(ibin,  hexpsig_z.GetBinContent(ibin + startbin-1))
        if hexpewk:
            hewk.SetBinContent(ibin,    hexpewk.GetBinContent(ibin + startbin-1))
        if hexpttbar:
            httbar.SetBinContent(ibin,  hexpttbar.GetBinContent(ibin + startbin-1))
        if hexpqcd:
            hqcd.SetBinContent(ibin,    hexpqcd.GetBinContent(ibin + startbin-1))

        hratio.SetBinContent(ibin, hexpfull.GetBinContent(ibin + startbin - 1))
        hratio.SetBinError(ibin,   hexpfull.GetBinError(ibin + startbin - 1))

        diff = horgdata.GetBinContent(ibin + startbin - 1) - hexpfull.GetBinContent(ibin + startbin - 1)
        # take the sigma as sqrt(data**2 + templates**2)
        # not 100% sure if this is the correct way to calculate pull
        sig = math.sqrt(horgdata.GetBinError(ibin + startbin - 1)**2 + hexpfull.GetBinError(ibin + startbin - 1)**2)
        hpull.SetBinContent(ibin, diff/(sig+1e-6))

    # deal with the uncertainty bar
    for ibin in range(1, hratio.GetNbinsX()+1):
        val = hratio.GetBinContent(ibin)
        err = hratio.GetBinError(ibin)
        hratio.SetBinContent(ibin, 1.0)
        if val!=0:
            hratio.SetBinError(ibin, err/val)
        else:
            hratio.SetBinError(ibin, 0.)

    hsig.SetFillColor(ROOT.TColor.GetColor(222, 90, 106))
    hsig.SetLineColor(1)

    hsig_z.SetFillColor(ROOT.TColor.GetColor(100, 192, 232))
    hsig_z.SetLineColor(1)

    hewk.SetFillColor(ROOT.TColor.GetColor("#E1F5A9"))
    hewk.SetLineColor(1)

    httbar.SetFillColor(ROOT.TColor.GetColor(155, 152, 204))
    httbar.SetLineColor(1)

    hqcd.SetFillColor(ROOT.TColor.GetColor(250, 202, 255))
    hqcd.SetLineColor(1)

    nevts = OrderedDict()
    nevts['data'] = hdata.Integral()
    nevts['sig'] = hsig.Integral()
    nevts['sig_z'] = hsig_z.Integral()
    nevts['ewk'] = hewk.Integral()
    nevts['ttbar'] = httbar.Integral()
    nevts['qcd'] = hqcd.Integral()

    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(1)
    hdata.SetMarkerColor(1)
    hdata.SetLineColor(1)
    hdata.Scale(1.0, "width")
    hqcd.Scale(1.0, "width")
    httbar.Scale(1.0, "width")
    hewk.Scale(1.0, "width")
    hsig.Scale(1.0, "width")
    hsig_z.Scale(1.0, "width")

    siglabels = {
        "muplus":  "W^{+}#rightarrow#mu^{+}#nu",
        "muminus": "W^{-}#rightarrow#mu^{-}#bar{#nu}",
        "mumu":    "Z#rightarrow #mu^{+}#mu^{-}"
        # "eplus":   "W^{+}#rightarrow e^{+}#nu",
        # "eminus":  "W^{-}#rightarrow e^{-}#bar{#nu}",
        # "ee":      "Z#rightarrow e^{+}e^{-}",
    }

    hs_gmc = ROOT.THStack("hs_stack_{}_{}".format(channel, suffix), "hs_stack")
    channel_forlabel = channel.replace("_pfmet", "").replace("_pfmt", "")
    labels_mc = []
    if nevts['qcd'] > 0.:
        hs_gmc.Add(hqcd)
        labels_mc.append("QCD")
    if nevts['ttbar'] > 0:
        hs_gmc.Add(httbar)
        labels_mc.append("t#bar{t}")
    if nevts['ewk'] > 0:
        hs_gmc.Add(hewk)
        labels_mc.append("EWK")
    if nevts['sig_z'] > 0:
        hs_gmc.Add(hsig_z)
        channel_z = channel_forlabel if ("plus" not in channel and "minus" not in channel) else channel_forlabel.replace("plus", "").replace("minus", "")*2
        labels_mc.append(siglabels[channel_z])
    if nevts['sig'] > 0:
        hs_gmc.Add(hsig)
        labels_mc.append(siglabels[channel_forlabel])
    labels_mc.reverse()

    if "ee" not in channel and "mumu" not in channel:
        # w's
        xlabel = None
        if "mt" in suffix:
            xlabel = "m_{T} [GeV]"
        elif "met" in suffix:
            xlabel = "MET [GeV]"
        outputname = f"{prepost}_w_{channel}_{suffix}"
    else:
        # z's
        xlabel = "m_{ll} [GeV]"
        outputname = f"{prepost}_z_{channel}_{suffix}"
    if x_label != "":
        xlabel = x_label

    ymaxs = {
        "muplus": 3.5e6,
        "muminus": 3.5e6,
        "mumu": 1.0e6,
    }

    ratiopanel_label = None
    if "_syst" in suffix:
        syst_name = suffix.split("_syst")[-1]
        if syst_name == "All":
            ratiopanel_label = "Total syst. unc."
        else:
            ratiopanel_label = f"Syst. unc. ({syst_name})"

    yrmin = 0.95
    yrmax = 1.05
    ypullmin = -3.99
    ypullmax = 3.99
    if "prefit" in prepost:
        yrmin = 0.85
        yrmax = 1.15
        ypullmin = -9.99
        ypullmax = 9.99

    drawconfigs = DrawConfig(
        xmin = bins.min(),
        xmax = bins.max(),
        xlabel = xlabel,
        ymin = 0,
        ymax = ymaxs[channel_forlabel] / bin_width,
        ylabel = f"Events / GeV",
        outputname = outputname,
        noCMS=False,
        dology=False,
        addOverflow=False,
        addUnderflow=False,
        yrmin=yrmin,
        yrmax=yrmax,
        yrlabel = "Obs/Pred"
    )
    DrawHistos(
        [
            hdata,
            hs_gmc
        ],
        ["Observed"]+labels_mc,
        drawconfigs.xmin,
        drawconfigs.xmax,
        drawconfigs.xlabel,
        drawconfigs.ymin,
        drawconfigs.ymax,
        drawconfigs.ylabel,
        drawconfigs.outputname,
        dology=drawconfigs.dology,
        dologx=drawconfigs.dologx,
        showratio=drawconfigs.showratio,
        yrmax = drawconfigs.yrmax,
        yrmin = drawconfigs.yrmin,
        yrlabel = drawconfigs.yrlabel,
        ypullmin=ypullmin,
        ypullmax=ypullmax,
        donormalize=drawconfigs.donormalize,
        ratiobase=drawconfigs.ratiobase,
        legendPos = drawconfigs.legendPos,
        redrawihist = drawconfigs.redrawihist,
        extraText = drawconfigs.extraText,
        noCMS = drawconfigs.noCMS,
        addOverflow = drawconfigs.addOverflow,
        addUnderflow = drawconfigs.addUnderflow,
        nMaxDigits = drawconfigs.nMaxDigits,
        hratiopanel=hratio,
        ratiopanel_label=ratiopanel_label,
        drawoptions=[
            'PE',
            'HIST same'
        ],
        showpull=showpull,
        hpulls=[hpull],
        W_ref = 600,
        is5TeV = False
    )

    ymaxs_logy = {
        "muplus": 3.5e9,
        "muminus": 3.5e9,
        "mumu": 1.0e9,
    }
    ymins_logy = {
        "muplus": 0.5e3,
        "muminus": 0.5e3,
        "mumu": 30,
    }

    drawconfigs = DrawConfig(
        xmin = bins.min(),
        xmax = bins.max(),
        xlabel = xlabel,
        ymin = ymins_logy[channel_forlabel],
        ymax = ymaxs_logy[channel_forlabel] / bin_width,
        ylabel = f"Events / GeV",
        outputname = outputname+"_log",
        dology=True,
        addOverflow=False,
        addUnderflow=False,
        yrmin=yrmin,
        yrmax=yrmax,
        yrlabel = "Obs/Pred"
    )
    DrawHistos(
        [
            hdata,
            hs_gmc
        ],
        ["Observed"]+labels_mc,
        drawconfigs.xmin,
        drawconfigs.xmax,
        drawconfigs.xlabel,
        drawconfigs.ymin,
        drawconfigs.ymax,
        drawconfigs.ylabel,
        drawconfigs.outputname,
        dology=drawconfigs.dology,
        dologx=drawconfigs.dologx,
        showratio=drawconfigs.showratio,
        yrmax = drawconfigs.yrmax,
        yrmin = drawconfigs.yrmin,
        yrlabel = drawconfigs.yrlabel,
        ypullmin=ypullmin,
        ypullmax=ypullmax,
        donormalize=drawconfigs.donormalize,
        ratiobase=drawconfigs.ratiobase,
        legendPos = drawconfigs.legendPos,
        redrawihist = drawconfigs.redrawihist,
        extraText = drawconfigs.extraText,
        noCMS = drawconfigs.noCMS,
        addOverflow = drawconfigs.addOverflow,
        addUnderflow = drawconfigs.addUnderflow,
        nMaxDigits = drawconfigs.nMaxDigits,
        hratiopanel=hratio,
        ratiopanel_label=ratiopanel_label,
        drawoptions=[
            'PE',
            'HIST same'
        ],
        showpull=showpull,
        hpulls=[hpull],
        W_ref = 600,
        is5TeV = False
    )

    return nevts


def GetPOIValue(ifilename, poiname = ""):
    """
    return the POI val and error given a postfit root file
    """
    f = ROOT.TFile(ifilename)
    tree = f.Get("fitresults")
    tree.GetEntry(0)
    val = getattr(tree, poiname)
    err = abs(getattr(tree, poiname+"_err"))
    return val, err


def ComparePOIs(vals_x: np.array, vals: list, errs: list, labels: list, colors: list, markers: list, output: str, is5TeV: bool):
    """
    compare the POI values with different selections
    """
    # print(vals_x)
    graphs = []
    nvals = len(vals)
    width = (vals_x[1]-vals_x[0])
    scale = 0.5
    for idx in range(nvals):
        val = vals[idx]
        err = errs[idx]
        color = colors[idx]
        marker = markers[idx]
        g = ROOT.TGraphErrors(len(vals_x), vals_x - scale*width/2. + idx*scale*width/(nvals-1.), val, np.zeros(len(vals_x)), err)
        g.SetLineColor(color)
        g.SetMarkerColor(color)
        g.SetMarkerStyle(markers[idx])
        graphs.append(g)
    ymin = 0.9
    ymax = 1.1
    w_var = None
    if "mt" in output:
        w_var = "m_{T} threshold [GeV]"
    elif "met" in output:
        w_var = "MET threshold [GeV]"
    DrawHistos(graphs, labels, vals_x[0]-width/2., vals_x[-1]+width/2., w_var, ymin, ymax, "POI / POI^{no cut}", output, dology=False, showratio=False, donormalize=False, drawoptions='EP', legendPos = [0.2, 0.7, 0.8, 0.8], noCMS = False, nMaxDigits = 3, legendNCols = 2, is5TeV = is5TeV, legendoptions=["LEP"]*nvals)


def result2json(ifilename: str, poiname: str, ofilename: str, hname: str = "nuisance_impact_mu"):
    """
    script to convert the postfit POI and impacts of nuisance parameters
    to json file, which will be used to make impact plots later
    """
    nameMap = {
        "Pol1shape": "QCD_pol1",
        "mcScale": "QCD_ScaledMC"
    }

    def getNuisName(nuis):
        result = nuis
        for key, val in nameMap.items():
            if nuis.endswith(key):
                #result = val
                result = nuis.replace(key, val)
                break
        if bool(re.match(r"\w*bin\d+shape", nuis)):
            result = ("QCD_" + nuis).replace("shape", "")
        return result.replace("lepEta_bin0_WpT_bin0_", "")

    ifile = ROOT.TFile(ifilename)
    himpact = ifile.Get(hname)
    tree = ifile.Get("fitresults")
    tree.GetEntry(0)

    # find the POI bin for poiname
    ibinX = -1
    for binX in range(1, himpact.GetNbinsX()+1):
        poi = himpact.GetXaxis().GetBinLabel(binX)
        if poi == poiname:
            ibinX = binX
            continue
    assert ibinX >=0, "Can not find the POI {} in the postfit file {}. Please check.".format(poiname, ifilename)

    results = OrderedDict()
    results['POIs'] = []
    val = getattr(tree, poiname)
    err = abs(getattr(tree, poiname+"_err"))
    poi = OrderedDict()
    poi['fit'] = [val-err, val, val+err]
    poi['name'] = poiname
    results['POIs'].append(poi)

    results['method'] = 'default'
    results['params'] = []

    # dump impacts
    impacts = OrderedDict()
    for ibinY in range(1, himpact.GetNbinsY()+1):
        nuis = himpact.GetYaxis().GetBinLabel(ibinY)
        impacts[nuis] = himpact.GetBinContent(ibinX, ibinY)

    # sort impacts, descending
    impacts = OrderedDict(sorted(list(impacts.items()), key=lambda x: abs(x[1]), reverse=True))

    pulls = OrderedDict()
    for nuis in list(impacts.keys()):
        val = getattr(tree, nuis)
        err = getattr(tree, nuis+"_err")
        err = abs(err)
        pulls[nuis] = [val - err, val, val + err]

    # save to results
    for nuis in list(impacts.keys()):
        systematic = OrderedDict()
        systematic['fit'] = pulls[nuis]
        systematic['groups'] = []
        systematic['impact_' + poiname] = impacts[nuis]
        systematic['name'] = getNuisName(nuis)
        systematic['prefit'] = [-1.0, 0., 1.0]
        systematic[poiname] = [poi['fit'][1] - impacts[nuis], poi['fit'][1], poi['fit'][1] + impacts[nuis]]
        systematic['type'] = "Gaussian"
        # print((getNuisName(nuis), pulls[nuis][1], pulls[nuis][1]-pulls[nuis][0], impacts[nuis]))

        results['params'].append(systematic)

    with open(ofilename, 'w') as fp:
        json.dump(results, fp, indent=2)


def DumpGroupImpacts(ifilename: str, poiname: str, hname = "nuisance_group_impact_mu"):
    """
    print out the grouped impacts
    """
    val_poi, err_poi = GetPOIValue(ifilename, poiname)

    ifile = ROOT.TFile(ifilename)
    himpact_grouped = ifile.Get(hname)

    # find the POI bin for poiname
    ibinX = -1
    for binX in range(1, himpact_grouped.GetNbinsX()+1):
        poi = himpact_grouped.GetXaxis().GetBinLabel(binX)
        if poi == poiname:
            ibinX = binX
            break
    assert ibinX >=0, "Can not find the POI {} in the postfit file {}. Please check.".format(poiname, ifilename)

    impacts = OrderedDict()
    for ibinY in range(1, himpact_grouped.GetNbinsY()+1):
        nuis = himpact_grouped.GetYaxis().GetBinLabel(ibinY)
        impacts[nuis] = himpact_grouped.GetBinContent(ibinX, ibinY) * 100.0 / val_poi

    stat_unc = impacts["stat"] * val_poi / 100.
    lumi_unc = 0.00 * val_poi if "ratio" not in poiname else 0.

    print("")
    print("#"*50)

    # adding BBB unc. to syst, not stats!
    err_poi = np.sqrt((impacts["binByBinStat"] / 100)**2 + (err_poi/val_poi)**2) * val_poi

    if lumi_unc > 0.:
        print(f"{ifilename:50s}|{poiname:30s}| poi = {val_poi:5.5f} +/- {stat_unc:5.5f} (stat) +/- {err_poi:5.5f} (syst) +/- {lumi_unc:5.5f} (lumi)")
    else:
        print(f"{ifilename:50s}|{poiname:30s}| poi = {val_poi:5.5f} +/- {stat_unc:5.5f} (stat) +/- {err_poi:5.5f} (syst)")

    # sort impacts, descending
    impacts = OrderedDict(sorted(list(impacts.items()), key=lambda x: abs(x[1]), reverse=True))

    print(f"\nPrint grouped nuisance impacts for {poiname} in {ifilename}")
    for nuis in list(impacts.keys()):
        print(f"{nuis:20}: {impacts[nuis]:.3f}")
    print()

    return impacts
