#!/usr/bin/env bash

export iteration='v14Asimov'
export fit_variable='pfmt_corr'
export datacard_name='card_mu'
export datacard_name_noxs="${datacard_name}_noxs"

export mass="120.0"

export gof_file_data="higgsCombineTest.GoodnessOfFit.mH120.root"
export gof_file_toys="higgsCombineTest.GoodnessOfFit.mH120.12345.root"


###########################
### choose ONE scenario ###
###########################

export gof_algo="saturated"

# export gof_algo="KS"

# export gof_algo="AD"

###########################
###########################
###########################


unset PYTHONHOME PYTHON_INCLUDE_PATH PYTHONIOENCODING PYTHONPATH PYTHON_VERSION
source combine_harvester_env.sh

cd cards/${iteration}/${fit_variable}/scan_wbin0_systAll/

combineTool.py -M CollectGoodnessOfFit --input ${gof_file_data} ${gof_file_toys} -m ${mass} -o gof.json
plotGof.py gof.json --statistic ${gof_algo} --mass ${mass} -o gof_plot --title-right=""

cd -
