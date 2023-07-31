#!/bin/bash

cd ./cards/v7_2/pfmt_corr/scan_wbin0_systAll


combineCards.py mumu=datacard_mumu.txt mumu_xsec=datacard_mumu_xsec_InAcc.txt > card_mu.txt


# syntax for absolute xsec, charge asymmetry, and xsec ratios



echo "Wplus_sig sumGroup = " >> card_mu.txt
echo "Wminus_sig sumGroup = " >> card_mu.txt
echo "Winc_sig sumGroup = " >> card_mu.txt
echo "Zinc_sig sumGroup = leplep_sig" >> card_mu.txt
echo "WchgAsym chargeMetaGroup = Wplus_sig Wminus_sig" >> card_mu.txt
echo "WchgRatio ratioMetaGroup = Wplus_sig Wminus_sig" >> card_mu.txt
echo "WplusZRatio ratioMetaGroup = Wplus_sig Zinc_sig" >> card_mu.txt
echo "WminusZRatio ratioMetaGroup = Wminus_sig Zinc_sig" >> card_mu.txt
echo "WZRatio ratioMetaGroup = Winc_sig Zinc_sig" >> card_mu.txt


text2hdf5.py card_mu.txt --maskedChan mumu_xsec --X-allow-no-background
combinetf.py card_mu.hdf5 --binByBinStat --computeHistErrors --saveHists --doImpacts --output card_mu.root --nThreads=24
cd -
