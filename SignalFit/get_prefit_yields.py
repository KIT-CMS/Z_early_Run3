import ROOT

filename = "root/earlyRun3_2022_combined.root"

histos = [
"VV#mm_corr-355862-357482-VV#Nominal#m_toFit",
"ST#mm_corr-355862-357482-ST#Nominal#m_toFit",
"DYtau#mm_corr-355862-357482-DYtau#Nominal#m_toFit",
"Wtau#mm_corr-355862-357482-Wtau#Nominal#m_toFit",
"DY#mm_corr-355862-357482-DY#Nominal#m_toFit",
"TT#mm_corr-355862-357482-TT#Nominal#m_toFit",
"data#mm_corr-355862-357482-data#Nominal#m_toFit",
"W#mm_corr-355862-357482-W#Nominal#m_toFit",
"VV#mmet_corr-355862-357482-VV#Nominal#pfmt_corr_neg",
"VV#mmet_corr-355862-357482-VV#Nominal#pfmt_corr_pos",
"ST#mmet_corr-355862-357482-ST#Nominal#pfmt_corr_neg",
"ST#mmet_corr-355862-357482-ST#Nominal#pfmt_corr_pos",
"DYtau#mmet_corr-355862-357482-DYtau#Nominal#pfmt_corr_neg",
"DYtau#mmet_corr-355862-357482-DYtau#Nominal#pfmt_corr_pos",
"Wtau#mmet_corr-355862-357482-Wtau#Nominal#pfmt_corr_neg",
"Wtau#mmet_corr-355862-357482-Wtau#Nominal#pfmt_corr_pos",
"DY#mmet_corr-355862-357482-DY#Nominal#pfmt_corr_neg",
"DY#mmet_corr-355862-357482-DY#Nominal#pfmt_corr_pos",
"TT#mmet_corr-355862-357482-TT#Nominal#pfmt_corr_neg",
"TT#mmet_corr-355862-357482-TT#Nominal#pfmt_corr_pos",
"data#mmet_corr-355862-357482-data#Nominal#pfmt_corr_neg",
"data#mmet_corr-355862-357482-data#Nominal#pfmt_corr_pos",
"W#mmet_corr-355862-357482-W#Nominal#pfmt_corr_neg",
"W#mmet_corr-355862-357482-W#Nominal#pfmt_corr_pos",
]

yield_dict = {}

ifile = ROOT.TFile(filename, "READ")

for i in histos:
    h_tmp_int = ifile.Get(i).Integral()
    print(f"{i}, {h_tmp_int}")
    yield_dict[i] = h_tmp_int

ifile.Close()

print()
print(" "*12+"Z"+" "*12+"W-"+" "*12+"W+")

histos.sort()

tmp_i = 0
for i in histos:

    if tmp_i == 0:
        print(i.split('#')[0], end="")
        print(" "*8, end="")


    print(yield_dict[i], end="")
    print(" "*8, end="")

    tmp_i = tmp_i + 1
    if tmp_i == 3:
        tmp_i = 0
        print()



val = {}
val['Z'] = {}
val['Wp'] = {}
val['Wm'] = {}

for i in histos:
    if i.endswith('toFit'):
        val['Z'][i.split('#')[0]] = yield_dict[i]
    if i.endswith('pos'):
        val['Wp'][i.split('#')[0]] = yield_dict[i]
    if i.endswith('neg'):
        val['Wm'][i.split('#')[0]] = yield_dict[i]

print()
print()
print()

print(val)


print()
print()
print()

print("\\begin{table}[htbp]")
print("\\topcaption{Expected yields from various processes in $\\PW^{+}$, $\\PW^{-}$, and $\\PZ$ muon final states. }")
print("\\centering")
print("\\begin{tabular}{lrrr}")
print("\\hline")
#print("Process &      $\\mathrm{W}^{+}\\rightarrow\\mu^{+}\\nu$ &     $\\mathrm{W}^{-}\\rightarrow\\mu^{-}\\bar{\\nu}$ &       $\\mathrm{Z}\\rightarrow\\mu^{+}\\mu^{-}$ \\\\")
print("Process &      $\\mu^{+}\\nu \, (\\mathrm{W}^{+})$ &     $\\mu^{-}\\bar{\\nu} \, (\\mathrm{W}^{-})$ &       $\\mu^{+}\\mu^{-} \, (\\mathrm{Z})$ \\\\")
print("\\hline")
print("W                                     &  " + '{:>12,.0f}'.format(val['Wp']['W']).replace(',','\\,') + " & " + '{:>12,.0f}'.format(val['Wm']['W']).replace(',','\\,') + "  & {} \\\\")
print("Z                                    &  " + '{:>12,.0f}'.format(val['Wp']['DY']).replace(',','\\,') + " & " + '{:>12,.0f}'.format(val['Wm']['DY']).replace(',','\\,') + "  & " + '{:>12,.0f}'.format(val['Z']['DY']).replace(',','\\,') + " \\\\")
print("EWK: & & & \\\\")
print("\\quad W ($\\rightarrow \\tau \\nu$)  &  " + '{:>12,.0f}'.format(val['Wp']['Wtau']).replace(',','\\,') + " & " + '{:>12,.0f}'.format(val['Wm']['Wtau']).replace(',','\\,') + "  & {} \\\\")
print("\\quad Z ($\\rightarrow \\tau \\tau$) &  " + '{:>12,.0f}'.format(val['Wp']['DYtau']).replace(',','\\,') + " & " + '{:>12,.0f}'.format(val['Wm']['DYtau']).replace(',','\\,') + "  & " + '{:>12,.0f}'.format(val['Z']['DYtau']).replace(',','\\,') + " \\\\")
print("\\quad VV                             &  " + '{:>12,.0f}'.format(val['Wp']['VV']).replace(',','\\,') + " & " + '{:>12,.0f}'.format(val['Wm']['VV']).replace(',','\\,') + "  & " + '{:>12,.0f}'.format(val['Z']['VV']).replace(',','\\,') + " \\\\")
print("\\quad single top                     &  " + '{:>12,.0f}'.format(val['Wp']['ST']).replace(',','\\,') + " & " + '{:>12,.0f}'.format(val['Wm']['ST']).replace(',','\\,') + "  & " + '{:>12,.0f}'.format(val['Z']['ST']).replace(',','\\,') + " \\\\")
print("$t\\bar{t}$                           &  " + '{:>12,.0f}'.format(val['Wp']['TT']).replace(',','\\,') + " & " + '{:>12,.0f}'.format(val['Wm']['TT']).replace(',','\\,') + "  & " + '{:>12,.0f}'.format(val['Z']['TT']).replace(',','\\,') + " \\\\")
print("QCD                                   &   xxx &   xxx &    {}   \\\\")
print("\\hline")
print("data                                  &   xxx &   xxx &    xxx   \\\\")
print("data                                  &  " + '{:>12,.0f}'.format(val['Wp']['data']).replace(',','\\,') + " & " + '{:>12,.0f}'.format(val['Wm']['data']).replace(',','\\,') + "  & " + '{:>12,.0f}'.format(val['Z']['data']).replace(',','\\,') + " \\\\")
print("\\hline")
print("\\end{tabular}")
print("\\label{tab:yields_muon_prefit}")
print("\\end{table}")
print()
