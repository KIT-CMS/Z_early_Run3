#!/usr/bin/env bash

export iteration='v14Asimov'
export fit_variable='pfmt_corr'
export datacard_name='card_mu'
export datacard_name_noxs="${datacard_name}_noxs"

export my_seed="12345"
export n_toys="100"

#export model_specific_options="--maskedChan muplus_xsec --maskedChan muminus_xsec --maskedChan mumu_xsec --X-allow-no-background"
export model_specific_options=""

###########################
### choose ONE scenario ###
###########################

export gof_algo="saturated"
export gof_toys="--toysFreq"

# export gof_algo="KS"
# export gof_toys=""

# export gof_algo="AD"
# export gof_toys=""

###########################
###########################
###########################

unset PYTHONHOME PYTHON_INCLUDE_PATH PYTHONIOENCODING PYTHONPATH PYTHON_VERSION
source combine_env.sh

cd cards/${iteration}/${fit_variable}/scan_wbin0_systAll/

combineCards.py muplus=datacard_muplus_pfmt_corr.txt muminus=datacard_muminus_pfmt_corr.txt mumu=datacard_mumu.txt > ${datacard_name_noxs}.txt

combine -M GoodnessOfFit ${datacard_name_noxs}.txt ${model_specific_options} --algo=${gof_algo} ${gof_toys}
combine -M GoodnessOfFit ${datacard_name_noxs}.txt ${model_specific_options} --algo=${gof_algo} ${gof_toys} -t ${n_toys} -s ${my_seed}

cd -
