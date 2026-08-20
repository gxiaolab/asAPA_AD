#! /usr/bin/bash
#$ -cwd
#$ -N homerprepplot
#$ -o log/job.$JOB_NAME.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_rt=1:00:00,h_data=1G
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_MOT=${DIR_WORK}/HOMER
DIR_INFO=${DIR_MOT}/STEP0_INFO
DIR_ONE=${DIR_MOT}/STEP1_GETSEQ
DIR_TWO=${DIR_MOT}/STEP2_HOMER
DIR_OUT=${DIR_MOT}/STEP3_PREPPLOT

####################

echo "--- $(date) --- Prep ---"

db="8mer"

dataset=${DIR_TWO}/asapa_events_homer_${db}.txt

dataset_name=$(basename $dataset .txt)

seq_len=7

N=$((seq_len+1))

python ./plot_homer_motifs_prep.py -i $dataset -n $N -o ${DIR_OUT}

echo "--- $(date) --- Done ---"
