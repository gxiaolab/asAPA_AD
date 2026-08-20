#! /usr/bin/bash
#$ -cwd
#$ -N homerplot
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
module load R

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_MOT=${DIR_WORK}/HOMER
DIR_INFO=${DIR_MOT}/STEP0_INFO
DIR_ONE=${DIR_MOT}/STEP1_GETSEQ
DIR_TWO=${DIR_MOT}/STEP2_HOMER
DIR_THREE=${DIR_MOT}/STEP3_PREPPLOT
DIR_OUT=${DIR_MOT}/STEP4_PLOTHOMER

####################

echo "--- $(date) --- Plot ---"

for dataset in ${DIR_THREE}/*;
do
	dataset_name=$(basename $dataset .txt)
	out_prefix=${DIR_OUT}/$dataset_name
	Rscript ./plot_homer_motifs.R ${dataset} $out_prefix
done

echo "--- $(date) --- Done ---"
