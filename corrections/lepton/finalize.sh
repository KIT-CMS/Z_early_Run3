export filedir='/ceph/moh/CROWN_samples/Run3V01/ntuples/2022/'
export execdir='/work/jdriesch/phd/Z_early_Run3/corrections/lepton/'

cd $filedir
for i in *
do
	cd $execdir
	python lepton_corrections.py -F mm -E -R 3 -V v2 --overwrite -I ${filedir}${i}/ --finalize
	python lepton_corrections.py -F mmet -E -R 3 -V v2 --overwrite -I ${filedir}${i}/ --finalize
	python lepton_corrections.py -F ee -E -R 3 -V v2 --overwrite -I ${filedir}${i}/ --finalize
        python lepton_corrections.py -F emet -E -R 3 -V v2 --overwrite -I ${filedir}${i}/ --finalize
done
