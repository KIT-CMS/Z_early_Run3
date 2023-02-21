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
    h_tmp_int = round(ifile.Get(i).Integral(),2)
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
