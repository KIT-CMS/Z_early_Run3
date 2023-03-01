#!/usr/bin/env bash


export iteration='v20'
export fit_variable='pfmt_corr'
export datacard_name='card_mu'
export script_dir="$PWD"

mkdir -p cards root plots logs pois {cards,root,plots}/${iteration}/${fit_variable}/

unset PYTHONHOME PYTHON_INCLUDE_PATH PYTHONPATH PYTHON_VERSION
source setenv.sh
python PrepareFits.py ${iteration}

echo "PrepareFits done"

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


unset PYTHON27PATH PYTHON3PATH PYTHON_VALGRIND_SUPP
source setenv.sh
python MakePostFitPlots.py ${iteration} &> logs/log-${iteration}.log

grep -r "poi =" logs/log-${iteration}.log > pois/res_pois_${iteration}.txt
python DrawDatacardSysts.py root/${iteration}/output_shapes_pfmt_corr_wbin0_systAll.root plots/syst_shapes_${iteration}/ &> logs/log-${iteration}_DrawDatacardSysts.log
