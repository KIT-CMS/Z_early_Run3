#!/usr/bin/env bash


export iteration='v13Asimov'
export fit_variable='pfmt_corr'
export datacard_name='card_mu'
export script_dir="$PWD"

mkdir -p {cards,root}/${iteration} plots

source setenv.sh
python PrepareFits.py ${iteration}

unset PYTHONHOME PYTHON_INCLUDE_PATH PYTHONIOENCODING PYTHONPATH PYTHON_VERSION
source combine_env.sh

for card_dir in cards/${iteration}/${fit_variable}/scan_wbin*; do
    cd $card_dir
    echo "starting $card_dir ..."
    (source ${datacard_name}.sh &> ${datacard_name}.log) &
    sleep 2
    cd $script_dir
done

echo "waiting for fits to finish ..."
wait

# tail -f  cards/${iteration}/${fit_variable}/scan_wbin0_syst*/${datacard_name}.log
# ls -lthr cards/${iteration}/*/scan_wbin*

unset PYTHON27PATH PYTHON3PATH PYTHON_VALGRIND_SUPP
source setenv.sh
python MakePostFitPlots.py ${iteration} &> log-${iteration}.log

grep -r "poi =" log-${iteration}.log > res_pois_${iteration}.txt
