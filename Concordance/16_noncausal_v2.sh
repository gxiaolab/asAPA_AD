#! /usr/bin/bash
#$ -cwd
#$ -N non-causal.v2
#$ -o logv2/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e logv2/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 1
#$ -t 1-6:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python/2.7.18

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_CONC=${DIR_WORK}/Concordance
DIR_ASARP=${DIR_WORK}/ASARP/asarp_final
DIR_ASARPOUT=${DIR_WORK}/ASARP/ASARP
events=${DIR_ASARP}/path/to/significant/asAPA_events.txt
DIR_ONE=${DIR_CONC}/STEP1_PREPROC
DIR_TWO=${DIR_CONC}/STEP2_RDD
DIR_FOUR=${DIR_CONC}/STEP4_SISCORES
DIR_FIVE=${DIR_CONC}/STEP5_PEAKS
DIR_SIX=${DIR_CONC}/STEP6_CAUSAL
DIR_OUT=${DIR_CONC}/STEP7_NONCAUSAL

###################

echo "--- $(date) --- Detect ---"

params=${DIR_OUT}/params/v2_params.txt
PARAMS=($(awk "NR==6" $params))

file=${PARAMS[0]}
indir=${PARAMS[1]} 
suff=${PARAMS[2]}
outf=${PARAMS[3]}

python2 ./detect.non-causal.v2.py -r $file -i $indir -s $suff -o $outf

echo "--- $(date) --- Done ---"
