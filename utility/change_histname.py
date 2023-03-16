import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--mtofit", action='store_true', default=False)
parser.add_argument("--corr_updn", action='store_true', default=False)

def corr_updn(f):

    tfile = ROOT.TFile(f, "UPDATE")

    processes = ["W", "ST", "VV", "DY", "TT", "Wtau", "DYtau", "data"]

    for p in processes:
        hn = f"{p}#mm_corr-355862-357482-{p}#Nominal#m_toFit_up;1"
        tfile.Delete(hn)
        hnup = f"{p}#mm_corr_up-355862-357482-{p}#Nominal#m_toFit_up"
        hup = tfile.Get(hnup)
        nup = hnup.replace("corr_up", "corr")
        hup.SetName(nup)
        hup.SetTitle(nup)
        hup.Write()


        hn = f"{p}#mm_corr-355862-357482-{p}#Nominal#m_toFit_dn;1"
        tfile.Delete(hn)
        hndn = f"{p}#mm_corr_dn-355862-357482-{p}#Nominal#m_toFit_dn"
        hdn = tfile.Get(hndn)
        ndn = hndn.replace("corr_dn", "corr")
        hdn.SetName(ndn)
        hdn.SetTitle(ndn)
        hdn.Write()
        
    for key in tfile.GetListOfKeys():
        name = key.GetName()
        if ("corr_up" in name) or ("corr_dn" in name):
            tfile.Delete(name+";1")
    tfile.Close()


def rename_mtofit(f):

    tfile = ROOT.TFile(f, "UPDATE")

    processes = ["W", "ST", "VV", "DY", "TT", "Wtau", "DYtau", "data"]

    for p in processes:
        hnup = f"{p}#mm_corr-355862-357482-{p}#Nominal#m_toFit_up"
        hup = tfile.Get(hnup)
        nup = hnup.replace("Nominal", "LepCorrUp").replace("m_toFit_up", "m_toFit")
        hup.SetName(nup)
        hup.SetTitle(nup)
        hup.Write()

        
        hndn = f"{p}#mm_corr-355862-357482-{p}#Nominal#m_toFit_dn"
        hdn = tfile.Get(hndn)
        ndn = hndn.replace("Nominal", "LepCorrDn").replace("m_toFit_dn", "m_toFit")
        hdn.SetName(ndn)
        hdn.SetTitle(ndn)
        hdn.Write()
    tfile.Close()


if __name__=="__main__":
    args = parser.parse_args()

    if args.mtofit:
        fin = "/work/jdriesch/earlyrun3/Z_early_Run3/output/earlyRun3_2022_combined.root"
        rename_mtofit(fin)

    if args.corr_updn:
        fin = "/work/jdriesch/earlyrun3/Z_early_Run3/output/earlyRun3_2022.root"
        corr_updn(fin)