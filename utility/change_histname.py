import ROOT

fin = "/work/jdriesch/earlyrun3/Z_early_Run3/output/earlyRun3_2022_combined.root"
fout = fin

tfile = ROOT.TFile(fin, "UPDATE")

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
