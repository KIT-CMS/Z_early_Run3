#!/usr/bin/env bash

export iteration='v14Asimov'
export fit_variable='pfmt_corr'
export datacard_name='card_mu'
export datacard_name_noxs="${datacard_name}_noxs"

unset PYTHONHOME PYTHON_INCLUDE_PATH PYTHONIOENCODING PYTHONPATH PYTHON_VERSION
source ../combine_env.sh

# create datacard for GoF tests

gof_main_dir="${PWD}"

cd ../cards/${iteration}/${fit_variable}/scan_wbin0_systAll/

combineCards.py muplus=datacard_muplus_pfmt_corr.txt muminus=datacard_muminus_pfmt_corr.txt mumu=datacard_mumu.txt > ${datacard_name_noxs}.txt

cp ${datacard_name_noxs}.txt ${gof_main_dir}/
cd $gof_main_dir



# prepare GoF tests

export gof_toy_dir="toys"
export gof_log_dir="logs"
mkdir -p $gof_toy_dir $gof_log_dir

# run data snapshot
echo "data..."

gof_algo="saturated"; gof_toys="--toysFreq"
./run_data.sh $gof_algo $gof_toys &> ${gof_log_dir}/data_${gof_algo}.log

# gof_algo="KS"; gof_toys=""
# ./run_data.sh $gof_algo $gof_toys &> ${gof_log_dir}/data_${gof_algo}.log

# gof_algo="AD"; gof_toys=""
# ./run_data.sh $gof_algo $gof_toys &> ${gof_log_dir}/data_${gof_algo}.log



# run toys
echo "toys..."

n_jobs=10
n_toys_per_job=10

start_seed=1
end_seed=`echo "$start_seed + $n_jobs -1" | bc`


for i_seed in $(seq $start_seed $end_seed); do

    gof_algo="saturated"; gof_toys="--toysFreq"
    ( ./run_toys.sh $gof_algo $gof_toys $n_toys_per_job $i_seed &> ${gof_log_dir}/toy_${gof_algo}_${n_toys_per_job}_${i_seed}.log ) &
    sleep 1s

    # gof_algo="KS"; gof_toys=""
    # ( ./run_toys.sh $gof_algo $gof_toys $n_toys_per_job $i_seed &> ${gof_log_dir}/toy_${gof_algo}_${n_toys_per_job}_${i_seed}.log ) &
    # sleep 1s

    # gof_algo="AD"; gof_toys=""
    # ( ./run_toys.sh $gof_algo $gof_toys $n_toys_per_job $i_seed &> ${gof_log_dir}/toy_${gof_algo}_${n_toys_per_job}_${i_seed}.log ) &
    # sleep 1s

done
