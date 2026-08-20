#! /usr/bin/bash
#$ -cwd
#$ -N peak_prep
#$ -o log2/job.$JOB_NAME.$TASK_ID.out
#$ -e log2/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_data=16G,h_rt=12:00:00
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
DIR_TWO=${DIR_CONC}/STEP2_RDD
DIR_FOUR=${DIR_CONC}/STEP4_SISCORES
DIR_OUT=${DIR_CONC}/STEP5_PEAKS

###################

echo "--- $(date) --- Start ---"

indir=${DIR_FOUR}/si.results
outdir=${DIR_OUT}/peaks.gmm
bin=10
mkdir -p $outdir
mkdir -p ${DIR_OUT}/files

logs=/path/to/log2/logs.txt

while read logfile
do
    log=$(echo -n "${logfile}" | cut -f1)
    grep "other gt" /path/to/log2/${log} |cut -d " " -f3|sort -u > ${DIR_OUT}/files/${log}_tmp || true
done < ${logs}

echo "--- $(date) --- Done ---"
