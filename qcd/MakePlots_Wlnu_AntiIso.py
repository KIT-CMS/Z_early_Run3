"""
test on the anti-Isolated region mixed with signal region
"""
import ROOT 
import numpy as np
from collections import OrderedDict
import sys, os , time
sys.path.append("/work/moh/tools/CMSPLOTS")
from myFunction import THStack2TH1

from SampleManager import DrawConfig, Sample, SampleManager

ROOT.TH1.AddDirectory(ROOT.kFALSE)
ROOT.TH2.AddDirectory(ROOT.kFALSE)
ROOT.TH3.AddDirectory(ROOT.kFALSE)
ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(32)

# boolean flag to set either muon or electron channel
# doMuon = False means the electron channel
doMuon = False

# boolean flag. if set to true, scale the MC cross section by 30%
applyScaling = False

def main():
    input_base = "inputs_EarlyRun3_V12"
    #input_base = "inputs_test"

    if not os.path.exists("root"):
        os.makedirs("root")

    print "Program start..."
    #time_init = time.time()

    # Test run
    #RUN = 323778
    #LUMI = 0.270935867e3

    #LUMI = 59.74e3
    LUMI = 59.832045316e3

    # Rochester correction
    ROOT.gInterpreter.ProcessLine('#include "RoccoR/RoccoR.cc"')
    ROOT.gInterpreter.ProcessLine('rc = RoccoR("RoccoR/RoccoR2018UL.txt")')
    ROOT.gRandom.SetSeed(1)

    # Trigger efficiency
    ROOT.gROOT.ProcessLine('TFile* f_trig_eff = TFile::Open("data/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root")')
    ROOT.gROOT.ProcessLine('TH1D* h_trig_eff_data  = (TH1D*)f_trig_eff->Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_charge_abseta_pt_efficiencyData")')
    ROOT.gROOT.ProcessLine('TH1D* h_trig_eff_mc  = (TH1D*)f_trig_eff->Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_charge_abseta_pt_efficiencyMC")')
    ROOT.gROOT.ProcessLine('f_trig_eff->Close();')

    if doMuon:
        #input_antiiso_data    = input_base+"/SingleMuon_Run2018D-UL2018/mmet/input_files.txt"
        input_antiiso_data    = input_base+"/SingleMuon_Run2018-UL2018/mmet/input_files.txt"
        input_antiiso_wl0     = input_base+"/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/mmet/input_files.txt"
        input_antiiso_wl1     = input_base+"/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/mmet/input_files.txt"
        input_antiiso_wl2     = input_base+"/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/mmet/input_files.txt"
        input_antiiso_zxx     = input_base+"/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/mmet/input_files.txt"
        input_antiiso_ttbar   = input_base+"/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X/mmet/input_files.txt"
        input_antiiso_tbar_1lep = input_base+"/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X/mmet/input_files.txt"

    else:
        input_antiiso_data    = input_base+"/EGamma_Run2018-UL2018/emet/input_files.txt"
        input_antiiso_wl0     = input_base+"/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/emet/input_files.txt"
        input_antiiso_wl1     = input_base+"/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/emet/input_files.txt"
        input_antiiso_wl2     = input_base+"/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/emet/input_files.txt"
        input_antiiso_zxx     = input_base+"/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X/emet/input_files.txt"
        input_antiiso_ttbar   = input_base+"/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X/emet/input_files.txt"
        input_antiiso_tbar_1lep = input_base+"/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X/emet/input_files.txt"

    ## for the QCD background estimation (data-driven)
    qcdnorm = 1.0
    mcscale = 1.0
    if applyScaling:
        # scale up the MC for 30%
        mcscale = 1.3
    DataAisoSamp  = Sample(input_antiiso_data, isMC=False, isMuon=doMuon, name="Data_aiso", isWSR=True, additionalnorm = qcdnorm, legend = 'QCD', color='226')
    
    # W -> lnu
    Wl0AisoSamp   = Sample(input_antiiso_wl0, isMC=True, isMuon=doMuon, name = "wl0_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale, lumi=LUMI)
    Wl1AisoSamp   = Sample(input_antiiso_wl1, isMC=True, isMuon=doMuon, name = "wl1_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale, lumi=LUMI)
    Wl2AisoSamp   = Sample(input_antiiso_wl2, isMC=True, isMuon=doMuon, name = "wl2_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale, lumi=LUMI)
    
    # ttbar
    TTbarAisoSamp  = Sample(input_antiiso_ttbar, isMC=True, isMuon=doMuon, name = "ttbar_dilepton_aiso",     isWSR=True, additionalnorm= qcdnorm * mcscale, lumi=LUMI)
    TT1LepAisoSamp = Sample(input_antiiso_tbar_1lep, isMC=True, isMuon=doMuon, name = "ttbar_1lepton_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale, lumi=LUMI)

    sampMan = SampleManager(DataAisoSamp, 
                                [Wl0AisoSamp, Wl1AisoSamp, Wl2AisoSamp, TTbarAisoSamp, 
                                TT1LepAisoSamp]
                            )
    sampMan.groupMCs(["ttbar_dilepton_aiso", "ttbar_1lepton_aiso"], "ttbar", 46, "t#bar{t}")
    sampMan.groupMCs(['wl0_aiso', 'wl1_aiso', 'wl2_aiso'], "wlnu", 92,"W#rightarrow#mu#nu")

    if doMuon:
        sampMan.DefineAll("Lep_pt",     "pt_rc_1")
    else:
        sampMan.DefineAll("Lep_pt",     "pt_1")
    sampMan.DefineAll("Lep_eta",    "eta_1")
    sampMan.DefineAll("Lep_phi",    "phi_1")
    sampMan.DefineAll("Lep_mass",   "mass_1")
    sampMan.DefineAll("Lep_q",      "q_1")

    # muon and electron isolation distributions are different
    # more coarse binning for electrons to make sure enough statistics
    if doMuon:
        #signal region
        sampMan.DefineAll("w_iso4", "(iso_1 < 0.15)")
        
        #background region
        sampMan.DefineAll("w_iso5", "(iso_1 > 0.20 && iso_1 < 0.25)")
        sampMan.DefineAll("w_iso6", "(iso_1 > 0.25 && iso_1 < 0.30)")
        sampMan.DefineAll("w_iso7", "(iso_1 > 0.30 && iso_1 < 0.35)")
        sampMan.DefineAll("w_iso8", "(iso_1 > 0.35 && iso_1 < 0.40)")
        sampMan.DefineAll("w_iso9", "(iso_1 > 0.40 && iso_1 < 0.45)")
        sampMan.DefineAll("w_iso10", "(iso_1 > 0.45 && iso_1 < 0.50)")
        sampMan.DefineAll("w_iso11", "(iso_1 > 0.50 && iso_1 < 0.55)")
        sampMan.DefineAll("w_iso12", "(iso_1 > 0.55 && iso_1 < 0.60)")
        sampMan.DefineAll("w_iso13", "(iso_1 > 0.60 && iso_1 < 0.65)")
        
        isobins = ["w_iso4", "w_iso5", "w_iso6", "w_iso7", "w_iso8", "w_iso9", "w_iso10", "w_iso11", "w_iso12", "w_iso13"]
    else:
        # separate the EB and EE and correct the electron isolation
        sampMan.DefineAll("eta_sc", "eta_1+deltaetaSC_1")
        sampMan.DefineAll("isEB",   "fabs(eta_sc) <= 1.479")
        
        #Signal Region
        #Signal region is defined differently for barrel and endcap
        sampMan.DefineAll("w_iso4", "isEB ? (iso_1 < 0.0478+0.506/pt_1) : (iso_1 < 0.0658+0.963/pt_1)")

        #Background Region
        sampMan.DefineAll("w_iso5", "(iso_1 > 0.20 && iso_1 < 0.25)")
        sampMan.DefineAll("w_iso6", "(iso_1 > 0.25 && iso_1 < 0.30)")
        sampMan.DefineAll("w_iso7", "(iso_1 > 0.30 && iso_1 < 0.35)")
        sampMan.DefineAll("w_iso8", "(iso_1 > 0.35 && iso_1 < 0.40)")
        sampMan.DefineAll("w_iso9", "(iso_1 > 0.40 && iso_1 < 0.45)")
        sampMan.DefineAll("w_iso10", "(iso_1 > 0.45 && iso_1 < 0.55)")
        sampMan.DefineAll("w_iso11", "(iso_1 > 0.55 && iso_1 < 0.70)")
        isobins = ["w_iso4", "w_iso5", "w_iso6", "w_iso7", "w_iso8", "w_iso9", "w_iso10", "w_iso11"]


    sampMan.DefineAll("met_phi", "metphi_uncorrected")
    sampMan.DefineAll("met", "met_uncorrected")

    #MMET:
    if doMuon:
        sampMan.DefineAll("mtCorr", "sqrt(2.*pt_rc_1*met_uncorrected*(1.-cos(phi_1 - metphi_uncorrected)))")
        sampMan.DefineAll("ptOmT", "pt_rc_1/mtCorr")
    else:
        sampMan.DefineAll("mtCorr", "sqrt(2.*pt_1*met_uncorrected*(1.-cos(phi_1 - metphi_uncorrected)))")
        sampMan.DefineAll("ptOmT", "pt_1/mtCorr")

    # charge
    lepname = "mu" if doMuon else "e"
    sampMan.DefineAll(lepname+"plus",  "Lep_q > 0")
    sampMan.DefineAll(lepname+"minus", "Lep_q < 0")
    chgbins = [lepname+"plus", lepname + "minus"]

    # WpT bins
    #sampMan.DefineAll("WpT_bin0",  "WpT>=0.")
    sampMan.DefineAll("WpT_bin0",  "Lep_pt>=0.")  #TODO: Changes this generic binning
    wptbins = ["WpT_bin0"]

    # eta bin: barral and endcap
    if doMuon:
        sampMan.DefineAll("lepEta_bin0", "1.0")
        sampMan.DefineAll("lepEta_bin1", "abs(Lep_eta) <= 1.2") 
        sampMan.DefineAll("lepEta_bin2", "abs(Lep_eta) > 1.2")
    else:
        sampMan.DefineAll("lepEta_bin0", "1.0")
        sampMan.DefineAll("lepEta_bin1", "abs(eta_sc) <= 1.44") 
        sampMan.DefineAll("lepEta_bin2", "abs(eta_sc) > 1.57")

    if doMuon:
        etabins = ["lepEta_bin0"]
        #etabins = ["lepEta_bin0", "lepEta_bin1", "lepEta_bin2"]
    else:
        etabins = ["lepEta_bin0"]
        #etabins = ["lepEta_bin0", "lepEta_bin1", "lepEta_bin2"]

    sampMan.DefineMC("Lep_trig_eff_data", "h_trig_eff_data->GetBinContent(\
        h_trig_eff_data->GetXaxis()->FindBin(Lep_q),\
        h_trig_eff_data->GetYaxis()->FindBin(TMath::Min((Float_t)2.3999, (Float_t)abs(Lep_eta))),\
        h_trig_eff_data->GetZaxis()->FindBin(TMath::Max((Float_t)26.001, TMath::Min((Float_t)199.999, (Float_t)Lep_pt)))\
    )")
    sampMan.DefineMC("Lep_trig_eff_mc", "h_trig_eff_mc->GetBinContent(\
        h_trig_eff_mc->GetXaxis()->FindBin(Lep_q),\
        h_trig_eff_mc->GetYaxis()->FindBin(TMath::Min((Float_t)2.3999, (Float_t)abs(Lep_eta))),\
        h_trig_eff_mc->GetZaxis()->FindBin(TMath::Max((Float_t)26.001, TMath::Min((Float_t)199.999, (Float_t)Lep_pt)))\
    )")
    sampMan.DefineMC("trig_eff_weight",
        "Lep_trig_eff_data / Lep_trig_eff_mc"
    )
    sampMan.DefineMC("weight_WoVpt", "weight_WoVpt_noTrigWeight*trig_eff_weight")

    for iso in isobins:
        for wpt in wptbins:
            for lepeta in etabins:
                for chg in chgbins:
                    sampMan.DefineAll("weight_{}_{}_{}_{}".format(chg, iso,  wpt, lepeta), "{} * weight_WoVpt * {} * {} * {}".format(iso, wpt, lepeta, chg))

    #print("sum of weight_WoVpt for data: ", sampMan.data.rdf.Sum("weight_WoVpt").GetValue())
    #print("Number of entries for data: ", sampMan.data.rdf.Count().GetValue())
    #return

    nbins = 24
    xmin = 0
    xmax = 120

    for iso in isobins:
        for wpt in wptbins:
            for lepeta in etabins:
                for chg in chgbins:
                    #print("at line 251, broadcasting", iso, wpt, lepeta, chg)
                    strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)

                    outputname = "histo_wjetsAntiIso_fullrange_mtcorr_" + strname 
                    sampMan.cacheDraw("mtCorr", outputname, nbins,  xmin, xmax, DrawConfig(xmin=xmin, xmax=xmax, xlabel="m_{T} [GeV]", dology=True, ymin=1,ymax=1e8, donormalizebin=False, addOverflow=False, addUnderflow=False, showratio=True), weightname = strname)
                    
                    outputname = "histo_wjetsAntiIso_fullrange_met_" + strname 
                    sampMan.cacheDraw("met", outputname, nbins,  xmin, xmax, DrawConfig(xmin=xmin, xmax=xmax, xlabel="MET [GeV]", dology=True, ymin=1,ymax=1e8, donormalizebin=False, addOverflow=False, addUnderflow=False, showratio=True), weightname = strname)

                    outputname = "histo_wjetsAntiIso_fullrange_ptOmT_" + strname 
                    sampMan.cacheDraw("ptOmT", outputname, nbins,  xmin, xmax, DrawConfig(xmin=xmin, xmax=xmax, xlabel="pT / mT [GeV]", dology=True, ymin=1,ymax=1e8, donormalizebin=False, addOverflow=False, addUnderflow=False, showratio=True), weightname = strname)

    strname = "weight_WoVpt"
    outputname = "histo_wjetsAntiIso_fullrange_lepEta_" + strname
    sampMan.cacheDraw("Lep_eta", outputname, 24, -2.4, 2.4, DrawConfig(xmin=-2.4, xmax=2.4, xlabel="Lepton Eta", dology=True, ymax=2e6, donormalizebin=False, addOverflow=False, addUnderflow=False, showratio=True), weightname = strname)

    print("Pre launch Draw")
    sampMan.launchDraw()
    print("Post launch Draw")

    if applyScaling:
        postfix = "nu_applyScaling.root"
    else:
        postfix = "nu.root"
    outfile = ROOT.TFile.Open("root/output_qcdshape_fullrange_"+lepname+"_"+postfix, "recreate")

    print("Computing Iso Mean")
    iso_mean = sampMan.data.rdf.Filter('iso_1<.15').Mean("iso_1").GetValue()
    print(240)
    iso_stat = ROOT.TStatistic("iso_mean")
    iso_stat.Fill(iso_mean, 1)
    #meta_tree = ROOT.TTree("meta", "meta")
    print(242)
    #meta_branch = meta_tree.Branch("iso_mean", 5)
    print(244)
    iso_stat.Write()
    #meta_tree.SetDirectory(outfile)
    print(246)
    #meta_tree.Write()
    #print("Iso mean = ", iso_mean.GetValue())
    print(249)

    for variable in ["mtcorr", "met", "ptOmT"]:
        hmts_comp = OrderedDict()
        for iso in isobins:
            for wpt in wptbins:
                for lepeta in etabins:
                    for chg in chgbins:
                        strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)

                        outputname = "histo_wjetsAntiIso_fullrange_"+variable+"_" + strname
                        hstacked = THStack2TH1(sampMan.hsmcs[outputname])
                        for ibin in xrange(hstacked.GetNbinsX()+1):
                            # hstacked should always be above 0
                            hstacked.SetBinContent(ibin, max(hstacked.GetBinContent(ibin), 0))
                        sampMan.hdatas[outputname].Add( hstacked, -1.0 )
                        hmts_comp[strname] = sampMan.hdatas[outputname]
                        hmts_comp[strname].SetName(outputname)

        odir = outfile.mkdir(variable)
        for wpt in wptbins:
            for iso in isobins[:-1]:
                for lepeta in etabins:
                    for chg in chgbins:
                        i = int(iso[5:])
                        iso_next = "w_iso" + str(i+1)
                        #strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)
                        #strname_next = "weight_{}_{}_{}_{}".format(chg, iso_next, wpt, lepeta)
                        strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)
                        strname_next = "weight_{}_{}_{}_{}".format(chg, iso_next, wpt, lepeta)
                        outputname = "histo_wjetsAntiIso_fullrange_"+ variable + "_" + strname

                        hcenter = hmts_comp[strname]
                        hup = hmts_comp[strname_next].Clone(outputname+"_shapeUp")
                        hdown = hmts_comp[strname_next].Clone(outputname+"_shapeDown")
                        hup.Scale(hcenter.Integral() / (hup.Integral()+1e-6))
                        for ibin in xrange(1, hcenter.GetNbinsX()+1):
                            center = hcenter.GetBinContent(ibin)
                            up = hup.GetBinContent(ibin)
                            hdown.SetBinContent(ibin, max(2*center - up, 0))

                            hcenter.SetBinContent(ibin, max(center, 0))
                            hup.SetBinContent(ibin, max(up, 0))
                        hcenter.SetDirectory(odir)
                        hcenter.Write()
                        hup.SetDirectory(odir)
                        hup.Write()
                        hdown.SetDirectory(odir)
                        hdown.Write()

    outfile.Close()


    sampMan.dumpCounts()

    print "Program end..."

    raw_input()
    
    return 

if __name__ == "__main__":
   main()
