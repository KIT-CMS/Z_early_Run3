#!/usr/bin/env bash

gof_algo="$1"
gof_toys="$2"

myname="_${gof_algo}"

#model_specific_options="--maskedChan muplus_xsec --maskedChan muminus_xsec --maskedChan mumu_xsec --X-allow-no-background"
model_specific_options=""

combine -M GoodnessOfFit -n $myname -m 125 -d ${gof_input} ${model_specific_options} --algo=${gof_algo} ${gof_toys}
mv higgsCombine${myname}.GoodnessOfFit.mH125.root ${gof_toy_dir}/
