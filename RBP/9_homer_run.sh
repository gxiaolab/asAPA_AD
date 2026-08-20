#! /usr/bin/bash
#$ -cwd
#$ -N homerrun
#$ -o log/job.$JOB_NAME.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_rt=2:00:00,h_data=4G
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
DIR_OUT=${DIR_MOT}/STEP2_HOMER

####################

echo "--- $(date) --- Homer ---"

seq_len=7
motif_size=$((seq_len+1))

db="8mer_preferred"

inf=${DIR_ONE}/higher_apa_events_${db}.fasta
bckg=${DIR_ONE}/lower_apa_events_${db}.fasta

homer2 denovo -i $inf -b $bckg -len $motif_size -strand + -o ${DIR_OUT}/asapa_events_homer_${db}.txt &>> ${DIR_OUT}/get_homer_motifs_${db}.out
homer2 denovo -i $bckg -b $inf -len $motif_size -strand + -o ${DIR_OUT}/asapa_events_homer_${db}_reverse.txt &>> ${DIR_OUT}/get_homer_motifs_${db}_reverse.out

echo "--- $(date) --- Done ---"
