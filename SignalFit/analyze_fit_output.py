import ROOT

fname = 'test_v20.root'
#fname = 'cards/v21/pfmt_corr/scan_wbin0_systAll/card_mu.root'

proceses = [
    'lepplus_sig',
    'lepminus_sig',
    'leplep_sig',
    'qcd_muplus',
    'qcd_muminus',
    'ewk',
    'tt',
]

fancy_process_names = {
    'lepplus_sig' : 'W+',
    'lepminus_sig': 'W-',
    'leplep_sig': 'Z',
    'qcd_muplus': 'QCD+',
    'qcd_muminus': 'QCD-',
    'ewk': 'EWK' ,
    'tt' : 'ttbar',
}


regions = ['muplus','muminus','mumu']
regionbins = {
    'muplus' : (1,20),
    'muminus' : (21,40),
    'mumu' : (41,70)
}


yields = {}
diff = {}

ifile = ROOT.TFile(fname, 'READ')
if ifile.IsOpen() == False:
    print("failed to open file ... abort")
    exit(1)

yields['data'] = {}
for t in ['prefit', 'postfit']:
    yields[t] = {}
    for r in regions:
        yields[t][r] = {}
        yields['data'][r] = {}
        yields['data'][r]['obs'] = ifile.Get('obs').Integral(regionbins[r][0], regionbins[r][1])
        for p in proceses:
            yields[t][r][p] = ifile.Get('expproc_'+p+'_'+t+'').Integral(regionbins[r][0], regionbins[r][1])


print(fname)
print()
print()

for r in regions:
    print(f'### {r} ###')
    print()
    diff[r] = {}
    print(f"Process           Prefit         Postfit      (Ratio)")
    print()
    for p in proceses:
        diff[r][p] = 1
        if yields['prefit'][r][p] > 0:
            diff[r][p] = yields['postfit'][r][p] / yields['prefit'][r][p]
            print(f"{fancy_process_names[p]:<10}: {yields['prefit'][r][p]:>12.2f} -> {yields['postfit'][r][p]:>12.2f}    ({diff[r][p]:>4.5f})")
        else:
            print(f"{fancy_process_names[p]:<10}: {yields['prefit'][r][p]:>12.2f} -> {yields['postfit'][r][p]:>12.2f}")
    print(f"{'data':<10}: {yields['data'][r]['obs']:>28.2f}")

    print()



for p in proceses:
    yields['prefit_'+p] = yields['prefit']['muplus'][p] + yields['prefit']['muminus'][p] +  yields['prefit']['mumu'][p]
    yields['postfit_'+p] = yields['postfit']['muplus'][p] + yields['postfit']['muminus'][p] +  yields['postfit']['mumu'][p]
    diff[p] = yields['postfit_'+p] / yields['prefit_'+p]

print('-'*20)
for p in proceses:
    print(f"{p:<14}: {yields['prefit_'+p]:>12.2f} -> {yields['postfit_'+p]:>12.2f} ({diff[p]:>5.8f})")
print()

mu = {}
itree = ifile.Get('fitresults')
itree.GetEntry(0)

mu['Wp'] = getattr(itree, 'lepplus_sig_mu')
mu['Wm'] = getattr(itree, 'lepminus_sig_mu')
mu['Z'] = getattr(itree, 'leplep_sig_mu')
mu['QCDp'] = getattr(itree, 'qcd_muplus_mu')
mu['QCDm'] = getattr(itree, 'qcd_muminus_mu')

print('-'*20)
for i in ['Wp', 'Wm', 'Z', 'QCDp', 'QCDm']:
    print(f"{i:<7}: {mu[i]:>43.8f}")
print()
