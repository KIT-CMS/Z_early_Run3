source utils/setup_cvmfs_sft.sh
source utils/setup_python.sh
export PYTHONPATH=$PWD:$PYTHONPATH

ERA=$1
INPUT=$2
VAR=$3
CHANNEL=$4
TAG=$5

python plotting/plot_shapes_control.py -l --era Run${ERA} --input ${INPUT} --variables ${VAR} --channels ${CHANNEL} --tag ${TAG}

