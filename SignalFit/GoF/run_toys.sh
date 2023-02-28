#!/usr/bin/env bash

gof_algo="$1"
n_toys="$2"
#my_seed="$3"
my_seed="-1"
gof_toys="$4"

echo $1 $2 $3 $4

myname="_${gof_algo}_${n_toys}_${my_seed}"

#model_specific_options="--maskedChan muplus_xsec --maskedChan muminus_xsec --maskedChan mumu_xsec --X-allow-no-background"
model_specific_options=""

combine -M GoodnessOfFit -n ${myname} -m 125 -d ${gof_input} ${model_specific_options} --algo=${gof_algo} ${gof_toys} -t ${n_toys} -s ${my_seed}
mv higgsCombine${myname}.GoodnessOfFit.mH125.${my_seed}.root ${gof_toy_dir}/
