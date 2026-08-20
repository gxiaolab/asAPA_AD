#! /usr/bin/bash
#$ -cwd
#$ -N peak_run
#$ -o log/job.$JOB_NAME.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_data=16G,h_rt=24:00:00
#$ -pe shared 4
#$ -t 1-1:1

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
DIR_OUT=${DIR_CONC}/STEP5_PEAKS

params=${DIR_OUT}/files/apa.peaks_r.job.params.txt
PARAMS=($(awk "NR==$SGE_TASK_ID" $params))

pref=${PARAMS[0]}
outf=${PARAMS[1]}
ref=/path/to/analysis/Concordance/STEP5_PEAKS/files/apa_r.other.gt.txt
minI=${PARAMS[2]}
nGMM=4
bin=${PARAMS[3]}
p=10

echo "--- $(date) --- Start ---"

python2 ./peak.si.rm.bg.downsample.py -i $pref -r $ref -m $minI -n $nGMM -o $outf -b $bin -p $p

echo "--- $(date) --- Done ---"
