source utils/setup_cvmfs_sft.sh
source utils/setup_python.sh
export PYTHONPATH=$PWD:$PYTHONPATH

ERA=$1
INPUT=$2
VAR=$3
CAT=$4
MATCH=$5
SEPVAR=$6
LABEL=$7
LATEX=$8
CHANNEL=$9
TAG=${10}

python plotting/plot_shapes_control_2018.py -l --era Run${ERA} --input ${INPUT} --variables ${VAR} --category-postfix ${CAT} --match-data ${MATCH} --seperate-variables ${SEPVAR} --lumi-label ${LABEL} --write-to-latex ${LATEX} --channels ${CHANNEL} --tag ${TAG}
