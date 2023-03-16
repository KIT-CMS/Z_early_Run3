import ROOT
import argparse
import numpy as np
from tqdm import tqdm
from utils import plot_ratio

parser = argparse.ArgumentParser()
parser.add_argument("--hists", action='store_true', default=False)
parser.add_argument("--weights", action='store_true', default=False)
parser.add_argument("--whists", action='store_true', default=False)
parser.add_argument("--plot", action='store_true', default=False)
parser.add_argument("--comp", action='store_true', default=False)

def makehists(fdict, rfile, hrange):
    rfile = ROOT.TFile(rfile, "RECREATE")
    for ch in fdict.keys():
        chain = ROOT.TChain("ntuple")
        for f in fdict[ch]:
            chain.Add(f)
            #print(f)
        rdf = ROOT.RDataFrame(chain)
        for var in hrange.keys():
            print(var)
            h = rdf.Histo1D(
                (ch+var, var, hrange[var][2], hrange[var][0], hrange[var][1]),
                var
            )
            print(h.Integral())
            h.Scale(1./h.Integral())
            h.Write()
    rfile.Close()
    return


def get_weights(rfile, hrange, rfile_w):
    hists = []
    for var in hrange.keys():
        # Get histograms from rootfile
        tf = ROOT.TFile(rfile, "READ")
        h_dt = tf.Get("data"+var)
        h_mc = tf.Get("mc"+var)
        h_dt.SetDirectory(ROOT.nullptr)
        h_mc.SetDirectory(ROOT.nullptr)
        tf.Close()

        h_weights = h_dt.Clone(var+"weight")

        for bin in range(hrange[var][2]+1):
            
            if h_dt.GetBinContent(bin+1)==0:
                weight = 0.
            elif h_mc.GetBinContent(bin+1)==0:
                weight = 1.
            else:
                weight = h_dt.GetBinContent(bin+1)/h_mc.GetBinContent(bin+1)
            
            h_weights.SetBinContent(bin+1, weight)

        hists.append(h_weights)

    tfw = ROOT.TFile(rfile_w, 'RECREATE')
    for h in hists:
        h.Write()
    tfw.Close()
    return
    

def make_weighted_hists(fdict, rfile, wfile, hrange):
    # Get weight hists from file
    wf = ROOT.TFile(wfile, "READ")
    wdict = {}
    for var in hrange.keys():
        wdict[str(var)] = wf.Get(var+"weight")
        wdict[var].SetDirectory(ROOT.nullptr)
    wf.Close()

    # Make file for weighted hists
    rf = ROOT.TFile(rfile, "RECREATE")

    # Go through different channels
    for ch in fdict.keys():

        # Open all root files from channel in rdf
        chain = ROOT.TChain("ntuple")
        for f in fdict[ch]:
            chain.Add(f)
        rdf = ROOT.RDataFrame(chain)
        rdf = rdf.Define('puweight', '0')
        quants = ["puweight"]
        
        # Go through different pileup quantities and change weights accordingly
        for var in hrange.keys():

            bins = np.linspace(hrange[var][0], hrange[var][1], hrange[var][2]+1)
            
            rdf = rdf.Define(var+'weight','1')

            if ch!='data':
                rdf = rdf.Redefine(var+'weight','0')
                for i in range(hrange[var][2]):
                    rdf = rdf.Redefine(
                        var+'weight',
                        'double w;'\
                        f'if ({var} >= {bins[i]} && {var} < {bins[i+1]})'\
                        f'w = {wdict[var].GetBinContent(i)};'\
                        f'else w = {var}weight;'\
                        'return w;'
                    )
                
            rdf = rdf.Redefine(
                'puweight', 
                f'puweight + 1./{len(hrange.keys())} * {var}weight'
                )
            quants += [var+'weight', var]

        # Check that the Integral remains the same after reweighting
        mean = rdf.Mean("puweight").GetValue()
        rdf = rdf.Redefine("puweight", f"puweight/{mean}")
        print(ch, rdf.Mean("puweight").GetValue())

        rdf.Snapshot("ntuple", f"{ch}.root", quants)
        
        for var in hrange.keys():
            print(var)
            h = rdf.Histo1D(
                (ch+var, var, hrange[var][2], hrange[var][0], hrange[var][1]),
                var,
                'puweight'
            )
            h.Scale(1./h.Integral())
            h.Write()

        if ch!='data':
            h_mm = rdf.Histo1D(
                ('pfmet_corr', "MET", 30, 0, 120),
                "pfmet_corr"
            )
            h_mm.Write()

            h_mm_w = rdf.Histo1D(
                ('pfmet_corr_weighted', "weighted MET", 30, 0, 120),
                "pfmet_corr",
                "puweight*"
            )
            h_mm_w.Write()


    rf.Close()
    return


def plot(rfile, hrange, outname):
    tf = ROOT.TFile(rfile, "READ")

    for var in hrange.keys():
        c = ROOT.TCanvas("c", "", 600, 600)
        #h_dt = tf.Get("data"+var)
        h_dt = tf.Get("data"+var)
        h_dy = tf.Get("mc"+var)
        
        plot_ratio(
            h_dt, 
            h_dy,
            labels = ['Data', 'MC'],
            title = var,
            outfile = f"plots/{var}{outname}.pdf",
            xrange = [hrange[var][0], hrange[var][1]],
            ratiorange = [.5, 1.5]
        )

    tf.Close()
    return


def make_hists_for_comparison(fdict):
    for ch in fdict.keys():
        chain = ROOT.TChain("ntuple")
        for f in fdict[ch]:
            chain.Add(f)
        rdf = ROOT.RDataFrame(chain)

        if ch == "wj":
            # uncorrected mt
            rdf = rdf.Define(
                "me",
                "sqrt(2*pt_1*pfmet_uncorrected*(1-cos(phi_1 - pfmetphi_uncorrected)))"
            )


            var = 'mt_1_uncorr'
            h = rdf.Histo1D(
                (var, "Transverse Mass", 20, 0, 120),
                var
            )
            h_w = rdf.Histo1D(
                (var+"_w", "Transverse Mass", 20, 0, 120),
                var,
                "puweight"
            )
            xrange = [0,120]
            title = "Transverse_Mass"
            rrange = [.9, 1.1]

        else:
            var = "m_vis"
            h = rdf.Histo1D(
                (var, "Di-Muon Mass", 30, 60, 120),
                var
            )
            h_w = rdf.Histo1D(
                (var+"_w", "Di-Muon Mass", 30, 60, 120),
                var,
                "puweight"
            )
            xrange = [60,120]
            title = "Visible_Mass"
            rrange = [.98, 1.02]

        tf = ROOT.TFile('tmp.root', "RECREATE")
        h.Write()
        h_w.Write()

        tf.Close()

        tf = ROOT.TFile('tmp.root', "READ")
        h = tf.Get(var)
        h_w = tf.Get(var+"_w")

        plot_ratio(
            h, 
            h_w,
            labels = ['MC', 'weighted MC'],
            title = title,
            outfile = f"plots/{title}.pdf",
            xrange = xrange,
            ratiorange = rrange      
        )
        tf.Close()



if __name__=='__main__':
    args = parser.parse_args()
    hists = 'hists/pt_hists.root'
    hists_w = 'hists/weighted_pt_hists.root'
    
    w_hists = 'weights.root'

    basedir = '/storage/9/jdriesch/earlyrun3/samples/Run3V06/ntuples_xsec_sf_scaleres_ptuncorr_pu_EraC/2022'
    #basedir =  '/storage/9/jdriesch/earlyrun3/samples/Run3V06/ntuples_xsec_sf_EraC/2022'
    muon = 'Muon_Run2022C-PromptReco-v1'
    #wj = 'WtoLNu_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X'
    dy = 'DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X'
    fdict = {
        'data': [
            f'{basedir}/Single{muon}/mm/Single{muon}_*.root',
            f'{basedir}/{muon}/mm/{muon}*.root',
        ],
        'mc': [
            f'{basedir}/{dy}/mm/{dy}*.root',
        ],
    }   
    """ 
    fdict_closure = {
        'wj': [
            f'{basedir}/{wj}/mmet/{wj}*.root'
        ]
    }
    """
    hrange = {
        'bosonpt': [0, 100, 100],
    }
    
    if args.hists:
        makehists(fdict, hists, hrange)

    if args.plot:
        plot(hists, hrange, "")

    if args.weights:
        get_weights(hists, hrange, w_hists)

    """
    if args.whists:
        make_weighted_hists(fdict, hists_w, w_hists, hrange)

    if args.comp:
        make_hists_for_comparison(fdict_closure)

    if args.plot:
        plot(hists_w, hrange, "_weighted")
    """