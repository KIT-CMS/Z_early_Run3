#!/usr/bin/env bash


###########################
### choose ONE scenario ###
###########################

export gof_algo="saturated"

# export gof_algo="KS"

# export gof_algo="AD"

###########################
###########################
###########################

export mass="125.0"

export gof_file_data="higgsCombine_${gof_algo}.GoodnessOfFit.mH125.root"
export gof_file_toys="higgsCombine_${gof_algo}_*_*.GoodnessOfFit.mH125.*.root"

# higgsCombine_saturated_100_1.GoodnessOfFit.mH125.1.root
# higgsCombine_saturated_100_2.GoodnessOfFit.mH125.2.root
# ...
# higgsCombine_saturated.GoodnessOfFit.mH125.root


# unset PYTHONHOME PYTHON_INCLUDE_PATH PYTHONIOENCODING PYTHONPATH PYTHON_VERSION
# source ../combine_harvester_env.sh

combineTool.py -M CollectGoodnessOfFit --input toys/${gof_file_data} toys/${gof_file_toys} -m ${mass} -o gof.json
plotGof.py gof.json --statistic ${gof_algo} --mass ${mass} -o gof_plot --title-right="" #--percentile 0 0.95
plotGof.py gof.json --statistic ${gof_algo} --mass ${mass} -o gof_plot_0p95 --title-right="" --percentile 0 0.95
plotGof.py gof.json --statistic ${gof_algo} --mass ${mass} -o gof_plot_0p90 --title-right="" --percentile 0 0.90

# cd -
