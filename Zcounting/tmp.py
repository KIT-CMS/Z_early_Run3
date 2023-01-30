lumi_online = 656.885772534
nZ_MC = 399510.87

nZ = 3.9752e+05
effHLT_var = 9.0057e-01
xsec_fid = 752.99
effIDIso = 0.92604156 * 0.99


lumi_Z = nZ * (lumi_online / nZ_MC)


lumi_Z = nZ / (xsec_fid * effIDIso*effIDIso)
xsec_Z = nZ / (lumi_online * effIDIso*effIDIso)

lumi_Z / lumi_online
xsec_Z / xsec_fid