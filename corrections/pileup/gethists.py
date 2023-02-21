import ROOT
import argparse
import json
import numpy as np
from tqdm import tqdm
from utils import plot_ratio

parser = argparse.ArgumentParser()
parser.add_argument("--hists", action='store_true', default=False)
parser.add_argument("--weights", action='store_true', default=False)
parser.add_argument("--whists", action='store_true', default=False)
parser.add_argument("--plot", action='store_true', default=False)

def makehists(fdict, rfile, hrange):
    rfile = ROOT.TFile(rfile, "RECREATE")
    for ch in fdict.keys():
        chain = ROOT.TChain("ntuple")
        for f in fdict[ch]:
            chain.Add(f)
        rdf = ROOT.RDataFrame(chain)
        for var in hrange.keys():
            print(var)
            h = rdf.Histo1D(
                (ch+var, var, hrange[var][2], hrange[var][0], hrange[var][1]),
                var
            )
            h.Scale(1./h.Integral())
            h.Write()
    rfile.Close()
    return

def get_weights(rfile, hrange):
    tf = ROOT.TFile(rfile, "READ")
    puweights = {}
    for var in hrange.keys():
        weights=[]
        h_dt = tf.Get("data"+var)
        h_dy = tf.Get("dy"+var)
        for bin in range(hrange[var][2]):
            if h_dy.GetBinContent(bin)==0:
                weight = 1.
            else:
                weight = h_dt.GetBinContent(bin)/h_dy.GetBinContent(bin)
            weights.append(weight)
        puweights[var] = weights
        # puweights[var] = [weights[i]/sum(weights)*len(weights) for i in range(len(weights))]
    with open('puweights.json', 'w') as f:
        json.dump(puweights, f)
    tf.Close()
    return
    

def make_weighted_hists(fdict, rfile, wfile, hrange):
    rfile = ROOT.TFile(rfile, "RECREATE")

    with open(wfile) as f:
        puweights = json.load(f)

    for ch in fdict.keys():
        chain = ROOT.TChain("ntuple")
        for f in fdict[ch]:
            chain.Add(f)
        rdf = ROOT.RDataFrame(chain)
        rdf = rdf.Define('puweight', '0')
        quants = ["m_vis", "puweight"]
        
        for var in hrange.keys():
            bins = np.linspace(hrange[var][0], hrange[var][1], hrange[var][2]+1)

            rdf = rdf.Define(var+'weight','1')
            
            if ch!='data':
                for i in range(hrange[var][2]):
                    rdf = rdf.Redefine(
                        var+'weight',
                        'double w;'\
                        f'if ({var} > {bins[i]} && {var} < {bins[i+1]})'\
                        f'w = {puweights[var][i]};'\
                        f'else w = {var}weight;'\
                        'return w;'
                    )
                
            rdf = rdf.Redefine(
                'puweight', 
                f'puweight + 1./{len(hrange.keys())} * {var}weight'
                )
            quants += [var+'weight', var]

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
                ('m_vis', "di-muon mass", 30, 60, 120),
                "m_vis"
            )
            h_mm.Write()

            h_mm_w = rdf.Histo1D(
                ('m_vis_weighted', "weighted di-muon mass", 30, 60, 120),
                "m_vis",
                "puweight"
            )
            h_mm_w.Write()


    rfile.Close()
    return


def plot(rfile, hrange, outname):
    tf = ROOT.TFile(rfile, "READ")

    for var in hrange.keys():
        c = ROOT.TCanvas("c", "", 600, 600)
        h_dt = tf.Get("data"+var)
        h_dy = tf.Get("dy"+var)
        
        plot_ratio(
            h_dt, 
            h_dy,
            labels = ['data', 'MC'],
            title = outname,
            outfile = f"plots/{var}{outname}.pdf",
            xrange = [hrange[var][0], hrange[var][1]],
            ratiorange = [0.5, 1.5]
        )

    if outname!="":
        c = ROOT.TCanvas("c", "", 600, 600)
        h_mm = tf.Get("m_vis")
        h_mm_w = tf.Get("m_vis_weighted")

        plot_ratio(
            h_mm, 
            h_mm_w,
            labels = ['MC', 'weighted MC'],
            title = outname,
            outfile = "plots/m_vis.pdf",
            xrange = [60,120],
            ratiorange = [0.98, 1.02]      
        )

    tf.Close()
    return


if __name__=='__main__':
    args = parser.parse_args()
    hists = 'hists/pileup_hists.root'
    whists = 'hists/weighted_pileup_hists.root'
    jsonfile = 'puweights.json'

    basedir = '/storage/9/jdriesch/earlyrun3/samples/Run3V06/ntuples/2022'
    muon = 'Muon_Run2022C-PromptReco-v1'
    dy = 'DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X'
    fdict = {
        'data': [
            f'{basedir}/Single{muon}/mm/Single{muon}_*.root',
            f'{basedir}/{muon}/mm/{muon}*.root',
        ],
        'dy': [
            f'{basedir}/{dy}/mm/{dy}*.root'
        ]
    }
    
    hrange = {
        'npvGood': [0, 60, 60],
        'rhoFastjetCentralChargedPileUp': [0, 50, 50],
        'rhoFastjetCentralCalo': [0, 35, 35]
    }
    
    if args.hists:
        makehists(fdict, hists, hrange)

    if args.plot:
        plot(hists, hrange, "")

    if args.weights:
        get_weights(hists, hrange)

    if args.whists:
        make_weighted_hists(fdict, whists, jsonfile, hrange)

    if args.plot:
        plot(whists, hrange, "_weighted")
