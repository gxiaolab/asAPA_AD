#! /usr/bin/bash
#$ -cwd
#$ -N rddtag
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=6:00:00,h_data=4G
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh

module load python

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_CONC=${DIR_WORK}/Concordance
DIR_ASARP=${DIR_WORK}/ASARP/asarp_final
DIR_ASARPOUT=${DIR_WORK}/ASARP/ASARP
events=${DIR_ASARP}/path/to/significant/asAPA_events.txt
DIR_ONE=${DIR_CONC}/STEP1_PREPROC
DIR_OUT=${DIR_CONC}/STEP2_RDD

snpfull=${DIR_ONE}/universal_tag_snp_rdd.bed
rnaseqids=${DIR_WORK}/ASARP/FQ/samples_names.txt

universal_tag=${DIR_ONE}/universal.tag.regionsfinal_clean.bed

####################

echo "--- $(date) --- Reformat ---"

rnaseqid=$(awk -v lin=${idx} 'FNR == (lin) {print $2;}' ${rnaseqids})

python3 ./tag.bed.rdd.py ${universal_tag} ${DIR_OUT}/${rnaseqid}.final.all_snvs.${rnaseqid}.rdd ${DIR_OUT}/${rnaseqid}.tag.rdd.bed

echo "--- $(date) --- Done ---"
