#!/usr/bin/env python

'''
create a plot for every systematic uncertainty in a datacard ROOT file with up, down and nominal distribution

usage: python DrawDatacardSysts.py datacard.root output
'''

#from __future__ import division
from ROOT import *
from sys import argv as cl
import re
import os
import sys

gStyle.SetOptStat(0)
gROOT.SetBatch(kTRUE)

region = []
process = []
syst = []
hlist = {}
n_plots = 0
n_warnings = 0

printint = False
debug = False

above = 0
below = 0

c = []
p1 = []
p2 = []

# create output folder if not present
if not os.path.isdir(cl[2]):
    os.mkdir(cl[2])

# open ROOT file and check if it exists
ifile = TFile(cl[1],"READ")
if ifile.IsOpen() == False:
    print("failed to open file ... abort")
    exit(1)
print('Creating plots from %s' %cl[1])

# create dictionary with histogram name as key and the histogram itself
for key in ifile.GetListOfKeys():
    if "PDF" in key.GetName():
        key.ReadObj().SetName(key.GetName())
    hlist[key.GetName()] = key.ReadObj()

# for key in ifile.GetListOfKeys():
#     hlist[key.ReadObj().GetName()] = key.ReadObj()

# loop over keys and figure out what regions, processes and systematics are present
for full_name in list(hlist):
    if full_name == '':
        continue
    n = '_'.join(full_name.split("_")[1:])
    tmp = n.split('_')[0]
    if not tmp in process and 'data_obs' not in tmp:
        process.append(tmp)
    tmp = '_'.join(n.split('_')[1:-1])
    if not tmp in region:
        region.append(tmp)
    tmp = n.split("_")[-1]
    tmp = re.sub('Up$', '', tmp)
    tmp = re.sub('Down$', '', tmp)
    if not tmp in syst and 'bin' not in tmp and 'nominal' not in tmp:
        syst.append(tmp)

# create the plots
for r in region:
    for p in process:
        for s in syst:

            sys.stdout.write('.')
            sys.stdout.flush()

            # figure out histogram names
            k_nom  = 'h_' + p + '_' + r + ''
            k_up   = 'h_' + p + '_' + r + '_' + s + 'Up'
            k_down = 'h_' + p + '_' + r + '_' + s + 'Down'

            h_name = 'h_' + p + '_' + r + '_' + s
            h_title = f'process: {p}, variable: {r}, systematic: {s}'

            # get histograms
            if k_up in list(hlist):
                h_nom  = hlist[k_nom]
                h_up   = hlist[k_up]
                h_down = hlist[k_down]
            else:
                continue
            # parameters for inconsistency checks
            yield_threshold = 2
            balance_threshold = 0.1
            reldiff_threshold = 0.5
            reldiff_upperlimit = 100
            ks_threshold = 1e-05

            # checks for yields
            printint = False
            if h_nom.Integral() > yield_threshold or h_up.Integral() > yield_threshold or h_down.Integral() > yield_threshold:

                balance = 1 - (((h_up.Integral() + h_down.Integral()) / 2) / h_nom.Integral())
                if abs(balance) > balance_threshold and h_nom.Integral() > yield_threshold*100:
                    #if not any(excl in s for excl in ('QCDscale','toppt','CMS_eff_')):
                    print('\n# Check\033[91m %s\033[39m -> unbalanced yield difference: %f' %(cl[2]+'/' + h_name + '.pdf', balance))
                    printint = True

                if h_up.Integral() == h_down.Integral():
                    if h_up.Integral() != h_nom.Integral():
                        if 'toppt' not in s:
                            print('\n# Check\033[91m %s\033[39m -> up and down variations are the same, but different compared to nominal:' %(cl[2]+'/' + h_name + '.pdf'))
                            printint = True

                    if h_up.Integral() == h_nom.Integral():
                        if not any(excl in s for excl in ('QCDscale','toppt','CMS_eff_')):
                            if not any(excl in p for excl in ('qcd')):
                                print('\n# Check\033[91m %s\033[39m -> up, down and nominal are the same:' %(cl[2]+'/' + h_name + '.pdf'))
                                printint = True
                if h_up.Integral()/h_nom.Integral()<0.99 and h_down.Integral()/h_nom.Integral()<0.99:
                    if 'toppt' not in s:
                        below += 1
                        print('\n# Check\033[91m %s\033[39m -> up, down are at least 1 percent below nominal' %(cl[2]+'/' + h_name + '.pdf'))
                        printint = True
                if h_up.Integral()/h_nom.Integral()>1.01 and h_down.Integral()/h_nom.Integral()>1.01:
                    if 'toppt' not in s:
                        above += 1
                        print('\n# Check\033[91m %s\033[39m -> up, down are at least 1 percent above nominal' %(cl[2]+'/' + h_name + '.pdf'))
                        printint = True

            if printint:
                print('# Integral nominal:\t\t %10.4f' %(h_nom.Integral()))
                print('# Integral up variation:\t %10.4f' %(h_up.Integral()))
                print('# Integral down variation:\t %10.4f' %(h_down.Integral()))
                n_warnings = n_warnings + 1


            # checks for distributions
            nbinsx = h_nom.GetNbinsX()
            found_bad_bins = False
            ks_up   = h_nom.KolmogorovTest(h_up)
            ks_down = h_nom.KolmogorovTest(h_down)
            bad_bins = []
            # check KS values as a starting point
            if ks_up < ks_threshold or ks_down < ks_threshold:

                for b in range(nbinsx):
                    # because ROOT
                    bin = b + 1
                    # check if bin is significant
                    if h_nom.GetBinContent(bin) > yield_threshold or h_up.GetBinContent(bin) > yield_threshold or h_down.GetBinContent(bin) > yield_threshold:

                        reldiff_up   = abs(1-h_up.GetBinContent(bin)/h_nom.GetBinContent(bin))
                        reldiff_down = abs(1-h_down.GetBinContent(bin)/h_nom.GetBinContent(bin))

                        # check if bin variation is within pre-defined region
                        if reldiff_upperlimit > reldiff_up > reldiff_threshold or reldiff_upperlimit > reldiff_down > reldiff_threshold:

                            diff_up   = abs(h_up.GetBinContent(bin)-h_nom.GetBinContent(bin))
                            diff_down = abs(h_down.GetBinContent(bin)-h_nom.GetBinContent(bin))

                            # check if covered by template uncertainty
                            if diff_up > (h_nom.GetBinError(bin)+h_up.GetBinError(bin)) or diff_down > (h_nom.GetBinError(bin)+h_down.GetBinError(bin)):
                                bad_bins.append(bin)
                                found_bad_bins = True

            bad_bins = [str(a) for a in bad_bins]
            if found_bad_bins:
                if not printint:
                    sys.stdout.write('\n')
                print('# Check\033[31m %s\033[39m -> bad bin(s):\033[93m %s \033[39m' %(cl[2]+'/' + h_name + '.pdf', ', '.join(bad_bins)))
                n_warnings = n_warnings + 1

            # normalize to max deviation
            maxbin = 0
            for h in [h_up,h_down,h_nom]:
                if h.GetMaximum() > maxbin:
                    maxbin = h.GetMaximum()
            h_up.SetMaximum( 1.2 * maxbin )

            # debug
            if debug:
                if any(fu in s for fu in ('scale_j','res_j','tH_UnclE', 'FSR_2017')):
                    print('%s ---> d:%f n:%f u:%f ---> diff_d:%f diff_u:%f' %(h_name, h_down.Integral(), h_nom.Integral(), h_up.Integral(), h_nom.Integral()-h_down.Integral(), h_nom.Integral()-h_up.Integral()))

            # canvas & style
            c.append(TCanvas(h_name,h_name,800,600))

            p1.append(TPad("pad1"+str(h_name)+str(n_plots), "pad1"+str(h_name)+str(n_plots), 0, 0.2, 1, 1.0))
            p1[-1].SetBottomMargin(0.05)
            p1[-1].SetRightMargin(0.05)
            p1[-1].SetLeftMargin(0.1)
            p1[-1].Draw()
            p1[-1].cd()


            h_nom.SetMinimum(0.001)
            h_up.SetMinimum(0.001)
            h_down.SetMinimum(0.001)

            h_nom.SetLineColor(1)
            h_up.SetLineColor(632)
            h_down.SetLineColor(600)

            h_up.SetTitle(h_title)

            h_up.GetYaxis().SetTitleSize(20)
            h_up.GetYaxis().SetTitleFont(43)
            h_up.GetYaxis().SetTitleOffset(1.55)

            h_up.GetYaxis().SetTitle("events")

            h_up.Draw()
            h_down.Draw('same')
            h_nom.Draw('same')

            # legend
            leg = TLegend(.71,.7,.95,.875)
            leg.SetBorderSize(0)
            leg.SetFillColor(0)
            leg.SetFillStyle(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.02)
            leg.SetTextAlign(12)
            int_nom = round(h_nom.Integral(),2)
            int_up = round(h_up.Integral(),2)
            int_down = round(h_down.Integral(),2)
            leg.AddEntry(h_nom,"nominal ("+str(int_nom)+")","L")
            leg.AddEntry(h_up,"up ("+str(int_up)+", "+str(round(int_up-int_nom,2))+")","L")
            leg.AddEntry(h_down,"down ("+str(int_down)+", "+str(round(int_down-int_nom,2))+")","L")
            leg.Draw()

            c[-1].cd()

            p2.append(TPad("pad2"+str(h_name), "pad2"+str(h_name), 0, 0.0, 1, 0.2))
            p2[-1].SetTopMargin(0.05)
            p2[-1].SetBottomMargin(0.05)
            p2[-1].SetRightMargin(0.05)
            p2[-1].SetLeftMargin(0.1)
            p2[-1].SetGridy()
            p2[-1].Draw()
            p2[-1].cd()

            h_diff_up = h_up.Clone(h_name+'diffup')
            h_diff_down = h_down.Clone(h_name+'diffdown')

            h_diff_up.Divide(h_nom)
            h_diff_down.Divide(h_nom)

            # h_diff_up.Sumw2()
            # h_diff_down.Sumw2()

            h_diff_up.SetStats(0)
            h_diff_down.SetStats(0)

            # ratio_max = max([h_diff_up.GetMaximum(5), h_diff_down.GetMaximum(5)])
            # print(ratio_max, ratio_min)

            ratio_max = 0
            for i in range(h_diff_up.GetNbinsX()):
                if abs(1 - h_diff_up.GetBinContent(i+1)) > ratio_max:
                    ratio_max = abs(1 - h_diff_up.GetBinContent(i+1))
                if abs(1 - h_diff_down.GetBinContent(i+1)) > ratio_max:
                    ratio_max = abs(1 - h_diff_down.GetBinContent(i+1))

            # ratio_max = abs(1 - h_diff_down.GetMaximum())

            h_diff_up.SetMaximum(1 + ratio_max*2)
            h_diff_up.SetMinimum(1 - ratio_max*2)
            # h_diff_up.SetMaximum(ratio_max)
            # h_diff_up.SetMinimum(ratio_min)

            h_diff_up.GetYaxis().SetNdivisions(505)
            h_diff_up.GetYaxis().SetLabelSize(0.125)
            h_diff_up.GetYaxis().SetTitle("ratio")

            h_diff_up.GetXaxis().SetLabelSize(0)

            h_diff_up.SetTitle('')

            h_diff_up.Draw('same')
            h_diff_down.Draw('same')


            # save and suppress write message
            gROOT.ProcessLine("gErrorIgnoreLevel = 2000;")
            c[-1].SaveAs(cl[2]+'/'+h_name+'.pdf')
            gROOT.ProcessLine("gErrorIgnoreLevel = -1;")
            n_plots = n_plots + 1


print()
print('DONE! (%i plots created), warnings: %i' %(n_plots, n_warnings))
print()
print('plots with up and down above nom: %i, below nom: %i' %(above,below))
print()
if below > above:
    print('more syst below nom than above! below: {}, above: {}'.format(below,above))
elif below < above:
    print('more syst above nom than below! below: {}, above: {}'.format(below,above))
print
