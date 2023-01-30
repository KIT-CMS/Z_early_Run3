#!/bin/bash

combineCards.py \
    muplus_pfmet=datacard_muplus_pfmet_corr.txt \
    muminus_pfmet=datacard_muminus_pfmet_corr.txt \
    muplus_pfmt=datacard_muplus_pfmt_corr.txt \
    muminus_pfmt=datacard_muminus_pfmt_corr.txt \
    mumu=datacard_mumu.txt \
    muplus_xsec=datacard_muplus_xsec_InAcc.txt \
    muminus_xsec=datacard_muminus_xsec_InAcc.txt \
    mumu_xsec=datacard_mumu_xsec_InAcc.txt \
    > card_mu.txt


# syntax for absolute xsec, charge asymmetry, and xsec ratios



echo "Wplus_sig sumGroup = lepplus_sig" >> card_mu.txt
echo "Wminus_sig sumGroup = lepminus_sig" >> card_mu.txt
echo "Winc_sig sumGroup = lepplus_sig lepminus_sig" >> card_mu.txt
echo "Zinc_sig sumGroup = leplep_sig" >> card_mu.txt
echo "WchgAsym chargeMetaGroup = Wplus_sig Wminus_sig" >> card_mu.txt
echo "WchgRatio ratioMetaGroup = Wplus_sig Wminus_sig" >> card_mu.txt
echo "WplusZRatio ratioMetaGroup = Wplus_sig Zinc_sig" >> card_mu.txt
echo "WminusZRatio ratioMetaGroup = Wminus_sig Zinc_sig" >> card_mu.txt
echo "WZRatio ratioMetaGroup = Winc_sig Zinc_sig" >> card_mu.txt


text2hdf5.py card_mu.txt --maskedChan muplus_xsec --maskedChan muminus_xsec --maskedChan mumu_xsec --X-allow-no-background
combinetf.py card_mu.hdf5 --binByBinStat --computeHistErrors --saveHists --doImpacts --output card_mu.root --nThreads=24
cd -