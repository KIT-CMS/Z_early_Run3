
####################

rm -rf root/test_v12Asimov/ cards/test_v12Asimov/* plots/*v12Asimov*

python PrepareFits.py v12Asimov

for card_dir in cards/test_v12Asimov/*/scan_wbin*; do eval "source $card_dir/card_mu.sh >&$card_dir/card_mu.log&"; done >&log&

tf cards/test_v12Asimov/pfmt_corr/scan_wbin0_syst*/card_mu.log
lt -r cards/test_v12Asimov/*/scan_wbin*

python MakePostFitPlots.py v12Asimov >&log-v12Asimov.log&

grep -r "poi =" log-v12Asimov.log >pois_v12Asimov.txt

####################

mkdir -p cards/test_v12Asimov/pfmet_pfmt
cd cards/test_v12Asimov/pfmet_pfmt
cp ../../../card_pfmet_pfmt.sh .

cp ../pfmet_corr/scan_wbin0_systAll/datacard_mu*pf* .
cp ../pfmet_corr/scan_wbin0_systAll/*xsec* .
cp ../pfmt_corr/scan_wbin0_systAll/datacard_mu*pf* .
cp ../pfmt_corr/scan_wbin0_systAll/datacard_mumu*.txt .

source card_pfmet_pfmt.sh >&card_mu.log&

####################

rsync -avzh --delete --progress moh@portal1.etp.kit.edu:/home/moh/CROWN_working/Z_early_Run3/SignalFit/plots .
rsync -avzh --delete --progress moh@portal1.etp.kit.edu:/home/moh/CROWN_working/Z_early_Run3/SignalFit/root .
rsync -avzh --delete --progress moh@portal1.etp.kit.edu:/home/moh/CROWN_working/Z_early_Run3/SignalFit/cards .

rsync -avzh --delete --progress /home/moh/CROWN_working/Z_early_Run3/SignalFit/root .
rsync -avzh --delete --progress /home/moh/CROWN_working/Z_early_Run3/SignalFit/cards .
