"""
code to generate the W(lnv) cards for tfCombine.
"""

from collections import OrderedDict
from re import A
import ROOT
import os
import json

class Process(object):
    """
    define one process and related information for the cards
    """
    def __init__(self,**kwargs):
        self.name = kwargs.get('name', 'sig')
        self.hname = kwargs.get('hname', 'htest')
        self.hsys  = kwargs.get('hsys', 'htest_sys')
        self.isObs = kwargs.get('isObs', False) # if it is the observation
        self.isSignal = kwargs.get('isSignal', False) # no contraint on the strength if true
        self.rate = kwargs.get('rate', -1.0)
        self.isMC = kwargs.get('isMC', False) # contraint on lumi if true
        self.isV = kwargs.get('isV', False) # recoil sys if true
        self.isQCD = kwargs.get('isQCD', False) # fake QCD
        self.xsecUnc = kwargs.get('xsecUnc', "0.10") # uncertainty on the cross section

        # a bit hacky, to point to the correct file location
        # - by default the fits run on lpc
        prefix_LPC = os.getcwd() + '/'
        self.fname = prefix_LPC + kwargs.get('fname', 'test.root')

        #  regulation
        if self.isObs:
            self.isSignal = False
            self.isMC = False
            self.isV = False
            self.isQCD = False
            self.xsecUnc = "-1.0"

        if self.isSignal:
            self.isMC = True # signal process is always uing MC template
            self.isV = True
            self.isQCD = False
            self.xsecUnc = "-1.0"

        if self.isQCD:
            self.isMC = False # QCD is always using data driven
            self.isV = False
            self.xsecUnc = "-1.0"


class Nuisance(object):
    """
    define one nuisance and related information for making cards
    """
    def __init__(self, **kwargs):
        self.name = kwargs.get('name', 'lumi')
        self.type = kwargs.get('type', 'lnN')
        assert self.type in ['lnN', 'shape'], "Nuisance type can only be lnN or shape"
        # valuemap saves how much one process is affected by the nuisance parameter
        # key is the process name, value is the uncertainty
        self.valuemap = kwargs.get('valuemap', {})

    def __getitem__(self, key):
        return self.valuemap.get(key, '-')

    def __setitem__(self, key, val):
        self.valuemap[key] = val


def WriteCard(data, processes, nuisgroups, cardname):
    """
    write the datacard provided with data, processes, and nuisances
    """
    dirpath = cardname.rpartition('/')[0]
    if not os.path.exists(dirpath):
        print(f"Make the directory {dirpath}")
        os.makedirs(dirpath)

    ofile = open(cardname, "w")
    ofile.write("imax 1 number of channels\n")
    ofile.write("jmax {} number of processes -1\n".format(len(processes)-1))
    ofile.write("kmax * number of nuisance parameters\n\n")
    # observation
    ofile.write("Observation -1\n\n")
    if data:
        # for signal mc xsec cards, no data entry
        ofile.write("shapes {pname:<30} * {fname:<40} {hname:<50}\n".format(pname = "data_obs", fname = data.fname, hname = data.hname))

    for proc in processes:
        ofile.write("shapes {pname:<30} * {fname:<40} {hname:<50} {hsys}$SYSTEMATIC\n".format(pname = proc.name, fname = proc.fname, hname = proc.hname, hsys = proc.hsys))
    ofile.write("\n")

    ofile.write("{:<40}".format("bin"))
    for proc in processes:
        # one data card is one bin
        ofile.write(" {binname:<15}".format(binname="bin1"))
    ofile.write("\n")

    ofile.write("{:<40}".format("process"))
    for proc in processes:
        ofile.write(" {pname:<15}".format(pname=proc.name))
    ofile.write("\n")

    ofile.write("{:<40}".format("process"))
    isig = 0
    ibkg = 1
    for proc in processes:
        if proc.isSignal or proc.isQCD:
            ofile.write(f" {isig:<15}")
            isig -= 1
        else:
            ofile.write(f" {ibkg:<15}")
            ibkg += 1
    ofile.write("\n")

    ofile.write("{:<40}".format("rate"))
    for proc in processes:
        ofile.write(" {:<15}".format("-1.0"))
    ofile.write("\n\n")

    ## write out all systematics
    for nuisgroup in nuisgroups.values():
        for nuis in nuisgroup:
            ofile.write(f"{nuis.name:<30} {nuis.type:<15}")
            for proc in processes:
                ofile.write(" {:<15}".format(nuis[proc.name]))
            ofile.write("\n")
    ofile.write("\n")

    # write out nuisance groups
    for gname, gnuisances in nuisgroups.items():
        ofile.write(f"{gname:<30} group =")
        for nuis in gnuisances:
            ofile.write(f" {nuis.name}")
        ofile.write("\n")

    ofile.close()


def MakeWJetsCards(fname_mc, fname_qcd, channel, rebinned = False, is5TeV = False, nMTBins = 9, outdir = "cards", applyLFU = False, fit_variable_w = "pfmet_corr", systs=None):
    # charge string for histogram names
    charge = "pos" if "plus" in channel else "neg"

    # prefix of all histo names
    prefix = ""
    if rebinned:
        prefix = "Rebinned_"

    # from lumi
    unc_lumi = {}
    unc_lumi['lumi_13p6TeV'] = 1.06

    # data
    data = Process(
        name = "data_obs", fname = fname_mc,
        hname = prefix + "h_data_{}_{}".format(fit_variable_w, charge),
        isObs = True,
    )

    # sig processes
    sigs = []
    sig = Process(
        name = GetSigName(channel, applyLFU), fname = fname_mc,
        hname = prefix + "h_w_{}_{}".format(fit_variable_w, charge),
        hsys  = prefix + "h_w_{}_{}_".format(fit_variable_w, charge),
        isSignal = True,
        isMC = True,
        isV = True,
        isQCD = False
    )
    sigs.append(sig)
    sig_z_channel = channel.replace("plus", "").replace("minus", "")*2
    sig_z = Process(
        name = GetSigName(sig_z_channel, applyLFU), fname = fname_mc,
        hname = prefix + "h_dy_{}_{}".format(fit_variable_w, charge),
        hsys  = prefix + "h_dy_{}_{}_".format(fit_variable_w, charge),
        isSignal = True,
        isMC = True,
        isV = True,
        isQCD = False
    )
    sigs.append(sig_z)

    # ttbar bkg
    ttbar = Process(
        name = "tt", fname = fname_mc,
        hname = prefix + "h_tt_{}_{}".format(fit_variable_w, charge),
        hsys  = prefix + "h_tt_{}_{}_".format(fit_variable_w, charge),
        isSignal = False,
        isMC = True,
        isV = False,
        isQCD = False,
        xsecUnc = "1.10"
    )

    # EWK bkg
    ewk = Process(
        name = "ewk", fname = fname_mc,
        hname = prefix + "h_ewk_{}_{}".format(fit_variable_w, charge),
        hsys  = prefix + "h_ewk_{}_{}_".format(fit_variable_w, charge),
        isSignal = False,
        isMC = True,
        isV = False,
        isQCD = False,
        xsecUnc = "1.10"
    )

    # QCD bkg
    qcd = Process(
        name = "qcd_"+channel, fname = fname_qcd,
        hname = prefix + "h_qcd_{}_{}".format(fit_variable_w, charge),
        hsys  = prefix + "h_qcd_{}_{}_".format(fit_variable_w, charge),
        isSignal = False,
        isMC = False,
        isV = False,
        isQCD = True,
    )

    # list of all processes
    processes = sigs + [ttbar, ewk, qcd]

    lepname = "mu" if "mu" in channel else "e"
    era = "13p6TeV"

    # define different nuisance parameters, their groups,
    # and the impacts on each process
    nuisgroups = OrderedDict()

    # nuis_lumi = Nuisance(name = "lumi_" + era, type = "lnN")
    # for proc in processes:
    #     if proc.isMC:
    #         nuis_lumi[proc.name] = unc_lumi[nuis_lumi.name]
    # nuisgroups["lumi"] = [nuis_lumi]

    nuisgroups["mcsec"] = []
    for proc in processes:
        if not proc.isSignal and proc.isMC:
            nuis_norm = Nuisance(name = "norm_" + proc.name, type = "lnN")
            nuis_norm[proc.name] = proc.xsecUnc
            nuisgroups["mcsec"].append(nuis_norm)

    # efficiency systematics
    nuisgroups["effsys"] = []
    for syst in ["SFTrk", "SFSta", "SFID", "SFIso", "SFTrg"]:
        nuisgroups["effsys"].append(Nuisance(name = syst, type = "shape"))
    for proc in processes:
        if not proc.isQCD:
            # all the samples except the QCD apply the corrections
            # so they should be affected by sf systematics
            for sysweight in nuisgroups["effsys"]:
                sysweight[proc.name] = 1.0

    # recoil correction systematics
    nuisgroups["recoil"] = []
    for syst in ["RecoilDoubleGauss", "RecoilSigOnlyFit", "RecoilZrap0", "RecoilZrap1", "RecoilZrap2"] + [f"RecoilStat{istat}" for istat in range(6)]:  # range(15)
        nuis_SysRecoil = Nuisance(name = syst, type = "shape")
        for proc in processes:
            if proc.isV:
                # all the V processes (including W and Z) apply the recoil corrections
                nuis_SysRecoil[proc.name] = 1.0
        nuisgroups["recoil"].append(nuis_SysRecoil)

    # qcd stat
    # this is hard coded for now. Will be improved later
    nuisgroups["qcdbkg"] = []
    nbins = nMTBins
    for syst in [f"bin{ibin}shape" for ibin in range(1, nbins+1)] + ["Pol1shape", "mcScale"]:
        nuis_QCDStat = Nuisance(name = syst, type = "shape")
        for proc in processes:
            if proc.isQCD:
                nuis_QCDStat[proc.name] = 1.0
        nuisgroups["qcdbkg"].append(nuis_QCDStat)

    # theory systematics
    nuisgroups["pdfscale"] = []
    for ipdf in range(1, 101):
        nuis_pdf = Nuisance(name = f"LHEPdfWeight{ipdf}", type = "shape")
        for proc in processes:
            if proc.isSignal:
                nuis_pdf[proc.name] = 1.0
        nuisgroups["pdfscale"].append(nuis_pdf)
    for syst in ["LHEPdfWeightAlphaS", "LHEScaleWeightMUF", "LHEScaleWeightMUR", "LHEScaleWeightMUFMUR"]:
        nuis_scale = Nuisance(name = syst, type = "shape")
        for proc in processes:
            if proc.isSignal:
                nuis_scale[proc.name] = 1.0
        nuisgroups["pdfscale"].append(nuis_scale)

    # clean up syst uncs for testing...
    nuisgroups_out = OrderedDict()
    if systs:
        assert type(systs) == list, systs
        for key in nuisgroups.keys():
            nuisgroups_out[key] = []
            for nuis in nuisgroups[key]:
                if nuis.name in systs:
                    nuisgroups_out[key].append(nuis)
            if len(nuisgroups_out[key]) < 1:
                del nuisgroups_out[key]
    else:
        nuisgroups_out = nuisgroups

    #
    # writing datacards
    #
    cardname = f"{outdir}/datacard_{channel}_{fit_variable_w}.txt"
    WriteCard(data, processes, nuisgroups_out, cardname)

    return cardname


def MakeZJetsCards(fname, channel, rebinned = False, is5TeV = False, outdir = "cards", applyLFU = False, systs=None):
    """
    Generate the combine datacard for Z+jets signal region
    """
    # prefix of all histo names
    prefix = ""
    if rebinned:
        prefix = "Rebinned_"

    # from lumi
    unc_lumi = {}
    unc_lumi['lumi_13p6TeV'] = 1.06

    # data
    data = Process(
        name = "data_obs", fname = fname,
        hname = prefix + "h_data_m_toFit",
        isObs = True,
    )

    # sig processes
    sigs = []
    sig_z = Process(
        name = GetSigName(channel, applyLFU), fname = fname,
        hname = prefix + "h_dy_m_toFit",
        hsys  = prefix + "h_dy_m_toFit_",
        isSignal = True,
        isMC = True,
        isV = True,
        isQCD = False
    )
    sigs.append(sig_z)

    # ttbar bkg
    ttbar = Process(
        name = "tt", fname = fname,
        hname = prefix + "h_tt_m_toFit",
        hsys  = prefix + "h_tt_m_toFit_",
        isSignal = False,
        isMC = True,
        isV = False,
        isQCD = False,
        xsecUnc = "1.10"
    )

    # EWK bkg
    ewk = Process(
        name = "ewk", fname = fname,
        hname = prefix + "h_ewk_m_toFit",
        hsys  = prefix + "h_ewk_m_toFit_",
        isSignal = False,
        isMC = True,
        isV = False,
        isQCD = False,
        xsecUnc = "1.10"
    )

    # list of all processes
    processes = sigs + [ttbar, ewk]

    lepname = "mu" if "mu" in channel else "e"
    era = "13p6TeV"

    # define different nuisance parameters, their groups,
    # and the impacts on each process
    nuisgroups = OrderedDict()

    # nuis_lumi = Nuisance(name = "lumi_" + era, type = "lnN")
    # for proc in processes:
    #     if proc.isMC:
    #         nuis_lumi[proc.name] = unc_lumi[nuis_lumi.name]
    # nuisgroups["lumi"] = [nuis_lumi]

    nuisgroups["mcsec"] = []
    for proc in processes:
        if not proc.isSignal and proc.isMC:
            nuis_norm = Nuisance(name = "norm_" + proc.name, type = "lnN")
            nuis_norm[proc.name] = proc.xsecUnc
            nuisgroups["mcsec"].append(nuis_norm)

    # efficiency systematics
    nuisgroups["effsys"] = []
    for syst in ["SFTrk", "SFSta", "SFID", "SFIso", "SFTrg"]:
        nuisgroups["effsys"].append(Nuisance(name = syst, type = "shape"))
    for proc in processes:
        if not proc.isQCD:
            # all the samples except the QCD apply the corrections
            # so they should be affected by sf systematics
            for sysweight in nuisgroups["effsys"]:
                sysweight[proc.name] = 1.0

    # theory systematics
    nuisgroups["pdfscale"] = []
    for ipdf in range(1, 101):
        nuis_pdf = Nuisance(name = f"LHEPdfWeight{ipdf}", type = "shape")
        for proc in processes:
            if proc.isSignal:
                nuis_pdf[proc.name] = 1.0
        nuisgroups["pdfscale"].append(nuis_pdf)
    for syst in ["LHEPdfWeightAlphaS", "LHEScaleWeightMUF", "LHEScaleWeightMUR", "LHEScaleWeightMUFMUR"]:
        nuis_scale = Nuisance(name = syst, type = "shape")
        for proc in processes:
            if proc.isSignal:
                nuis_scale[proc.name] = 1.0
        nuisgroups["pdfscale"].append(nuis_scale)

    # clean up syst uncs for testing...
    nuisgroups_out = OrderedDict()
    if systs:
        assert type(systs) == list, systs
        for key in nuisgroups.keys():
            nuisgroups_out[key] = []
            for nuis in nuisgroups[key]:
                if nuis.name in systs:
                    nuisgroups_out[key].append(nuis)
            if len(nuisgroups_out[key]) < 1:
                del nuisgroups_out[key]
    else:
        nuisgroups_out = nuisgroups

    #
    # writing datacards
    #
    cardname = f"{outdir}/datacard_{channel}.txt"
    WriteCard(data, processes, nuisgroups_out, cardname)

    return cardname


def combineCards(labels, cards, oname):
    """
    genrate and run the command to combine cards in different bins
    """
    cmd = "cd cards; combineCards.py"
    for label, card in zip(labels, cards):
        cmd += " {}={}".format(label, card)
    cmd += " > {}; cd ../".format(oname)
    print(("combine cards with {}".format(cmd)))
    os.system(cmd)


def GenerateRunCommand(output: str, cards: list, channels: list, cards_xsec: list = [], prefix: str = "./", applyLFU: bool = False, doAsimov: bool = False):
    """
    generate the script with commands to run combine.
    inputs can probably be improved here
    """
    assert len(cards) == len(channels), "cards and channels should have the same length"
    if len(cards_xsec) > 0:
        # if provided xsec cards, then it should have the same length as the hist cards and channels
        assert len(cards_xsec) == len(channels), "xsec cards and channels should have the same length"

    # these partitions can be changed depending on the directories and organizations
    card0 = cards[0]
    workdir = prefix + card0.rpartition('/')[0]

    for idx, card in enumerate(cards):
        cards[idx] = card.split('/')[-1]

    for idx, card_xsec in enumerate(cards_xsec):
        cards_xsec[idx] = card_xsec.split('/')[-1]

    cmd = ""
    cmd += "#!/bin/bash\n\n"
    # cmd += f"cd {workdir}\n\n\n"

    cmd += f"combineCards.py"
    for channel, card in zip(channels, cards):
        if card == None:
            continue
        cmd += f" {channel}={card}"
    if len(cards_xsec) > 0:
        for channel, card_xsec in zip(channels, cards_xsec):
            if card_xsec == None:
                continue
            cmd += f" {channel}_xsec={card_xsec}"
    cmd += f" > {output}.txt\n"

    if len(cards_xsec) > 0:
        wplus = set()
        wminus = set()
        zinc = set()
        # add the syntax to do the absolute xsec, charge asymmetry, and ratio measurements
        for channel in channels:
            signame = GetSigName(channel, applyLFU)
            if "plus" in channel:
                wplus.add(signame)
            elif "minus" in channel:
                wminus.add(signame)
            else:
                # others are z's
                zinc.add(signame)
        winc = wplus.union(wminus)

        cmd += '\n\n'
        cmd += '# syntax for absolute xsec, charge asymmetry, and xsec ratios\n'
        cmd += '\n\n\n'
        cmd += f'echo \"Wplus_sig sumGroup = {" ".join(sig for sig in wplus)}\" >> {output}.txt\n'
        cmd += f'echo \"Wminus_sig sumGroup = {" ".join(sig for sig in wminus)}\" >> {output}.txt\n'
        cmd += f'echo \"Winc_sig sumGroup = {" ".join(sig for sig in winc)}\" >> {output}.txt\n'
        cmd += f'echo \"Zinc_sig sumGroup = {" ".join(sig for sig in zinc)}\" >> {output}.txt\n'
        cmd += f'echo \"WchgAsym chargeMetaGroup = Wplus_sig Wminus_sig\" >> {output}.txt\n'
        cmd += f'echo \"WchgRatio ratioMetaGroup = Wplus_sig Wminus_sig\" >> {output}.txt\n'
        cmd += f'echo \"WplusZRatio ratioMetaGroup = Wplus_sig Zinc_sig\" >> {output}.txt\n'
        cmd += f'echo \"WminusZRatio ratioMetaGroup = Wminus_sig Zinc_sig\" >> {output}.txt\n'
        cmd += f'echo \"WZRatio ratioMetaGroup = Winc_sig Zinc_sig\" >> {output}.txt\n'
        cmd += '\n\n'

    cmd += f"text2hdf5.py {output}.txt"
    if len(cards_xsec) > 0:
        for channel in channels:
            cmd += f" --maskedChan {channel}_xsec"
        cmd += " --X-allow-no-background"
    cmd += "\n"
    if doAsimov:
        cmd += f"combinetf.py {output}.hdf5 --binByBinStat --computeHistErrors --saveHists --doImpacts --output {output}.root --nThreads=24 -t -1\n"
    else:
        cmd += f"combinetf.py {output}.hdf5 --binByBinStat --computeHistErrors --saveHists --doImpacts --output {output}.root --nThreads=24\n"
    # cmd += "cd -\n"

    # write outputs to scripts
    with open(f"{workdir}/{output}.sh", "w") as f:
        f.write(cmd)
    os.system(f"chmod +x {workdir}/{output}.sh")


def GetXSec(channel: str, is5TeV: bool = False):
    """
    return the MC Xsec given a channel and sqrtS.
    for longer term can read this from some root/txt files
    """

    f = open('../acceptance/acceptance.json', 'r')
    info = json.load(f)
    acc = "genlep_acc"  # "genlepPreFSR_acc"

    n_Wpos = info["Wpos"]["mu"]["genlepPreFSR_den_pw"] - info["Wpos"]["mu"]["genlepPreFSR_den_nw"]
    n_Wneg = info["Wneg"]["mu"]["genlepPreFSR_den_pw"] - info["Wneg"]["mu"]["genlepPreFSR_den_nw"]

    xsecs = {}
    xsecs["13p6TeV"] = {}
    xsecs["13p6TeV"]["muplus"] =  0.5 * 42133.267 * (n_Wpos/(n_Wpos+n_Wneg)) * info["Wpos"]["mu"][acc]
    xsecs["13p6TeV"]["muminus"] = 0.5 * 42133.267 * (n_Wneg/(n_Wpos+n_Wneg)) * info["Wneg"]["mu"][acc]
    xsecs["13p6TeV"]["mumu"] =    0.5 * 4147.5333 * info["Z"]["mu"][acc]

    xsecs["13p6TeV"]["eplus"] = 1e6
    xsecs["13p6TeV"]["eminus"] = 1e6
    xsecs["13p6TeV"]["ee"] = 1e6

    return xsecs["13p6TeV"][channel]


def MakeXSecCard(channel: str, is5TeV: bool = False, outdir: str = "cards", applyLFU: bool = False):
    """
    Generate the root file with xsec result, and also the datacard
    """
    # get the xsec result
    xsec = GetXSec(channel, is5TeV)

    # write the xsec result to a root file
    fname = f"{outdir}/xsec_{channel}.root"
    f = ROOT.TFile(fname, "RECREATE")
    hname = f"xsec_{channel}_inAcc"
    h = ROOT.TH1D(hname, hname, 1, 0, 1)
    h.SetBinContent(1, xsec)
    h.SetBinError(1, 0)
    h.Write()
    f.Close()

    # sig mc xsec
    sig = Process(name = GetSigName(channel, applyLFU), fname = fname,
                 hname = hname,
                 hsys = hname+"_",
                 isSignal = True,
                 isMC = True,
                 isV = True,
                 isQCD = False
                 )

    # get the datacard
    cardname = f"{outdir}/datacard_{channel}_xsec_InAcc.txt"

    nuisgroups = OrderedDict()
    processes = [sig]
    WriteCard(None, processes, nuisgroups, cardname)

    return cardname


def GetSigName(channel: str, applyLFU: bool = False):
    """
    return the signal process name
    """
    signame = channel
    if applyLFU:
        signame = signame.replace("e", "lep").replace("mu", "lep")
    return signame + "_sig"
