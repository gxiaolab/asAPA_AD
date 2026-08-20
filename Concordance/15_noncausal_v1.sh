#! /usr/bin/bash
#$ -cwd
#$ -N non-causal.v1
#$ -o logv1/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e logv1/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 1
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
DIR_FIVE=${DIR_CONC}/STEP5_PEAKS
DIR_SIX=${DIR_CONC}/STEP6_CAUSAL
DIR_OUT=${DIR_CONC}/STEP7_NONCAUSAL/v1

###################

echo "--- $(date) --- Set variables ---"

mkdir -p ${DIR_OUT}/AD
mkdir -p ${DIR_OUT}/Control

ref_dir=${DIR_SIX}/v1

ref_Control=${ref_dir}/Control
ref_AD=${ref_dir}/AD

outf_Control=${DIR_OUT}/Control
outf_AD=${DIR_OUT}/AD
mkdir -p ${outf_Control}
mkdir -p ${outf_AD}

suff_AD=AD/
suff_Control=Control/

echo "--- $(date) --- Get Non Causal ---"

for file in $ref_Control/*
do
        name=$(basename $file)
        tissue="Tissue Type"
        indir=${DIR_FOUR}/si.results/$tissue
        python2 ./detect.non-causal.py -r $file -i $indir -s $suff_Control -o $outf_Control/$name
done

for file in $ref_AD/*
do
        name=$(basename $file)
        tissue="Tissue Type"
        indir=${DIR_FOUR}/si.results/$tissue
        python2 ./detect.non-causal.py -r $file -i $indir -s $suff_AD -o $outf_AD/$name
done

echo "--- $(date) --- Done ---"

