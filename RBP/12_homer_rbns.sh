#! /usr/bin/bash
#$ -cwd
#$ -N motifScan
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=6:00:00,h_data=4G,highp
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
DIR_FOUR=${DIR_MOT}/STEP4_PLOTHOMER
DIR_OUT=${DIR_MOT}/STEP5_RBNS

####################

echo "--- $(date) --- Plot ---"

rbns=/path/to/rbns_data.txt
homer_motif_dir=${DIR_THREE}
outdir=${DIR_OUT}

indir=${homer_motif_dir}
outf=$DIR_OUT/events_rbns_overlap.txt

python ./scan_motifs_against_rbns.py -i $indir -r $rbns -o $outf

echo "--- $(date) --- Done ---"
