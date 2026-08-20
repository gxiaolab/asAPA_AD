#! /usr/bin/bash
#$ -cwd
#$ -N non-causal.param
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
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
DIR_OUT=${DIR_CONC}/STEP7_NONCAUSAL

###################

echo "--- $(date) --- Setup ---"

tmp_v2=${DIR_OUT}/v2_tmp.txt
tmp_v2b=${DIR_OUT}/v2b_tmp.txt

ls ${DIR_SIX}/v2/*/* >> ${tmp_v2}
ls ${DIR_SIX}/v2b/*/* >> ${tmp_v2b}

mkdir -p ${DIR_OUT}/v2/AD
mkdir -p ${DIR_OUT}/v2/Control

mkdir -p ${DIR_OUT}/v2b/AD
mkdir -p ${DIR_OUT}/v2b/Control

mkdir -p ${DIR_OUT}/params

echo "--- $(date) --- Create file: v2 ---"

while read file
do
        filename=$(echo -n "${file}" | cut -f 1)
        base=$(basename "$filename")
        tissue="DLPFC"

        ## Check Status ##
        if [[ "$filename" == *"Control"* ]]; then
                statuss="Control"
        else
                statuss="AD"
        fi
        indir=${DIR_FOUR}/si.results/${tissue}
        suffix="${statuss}/"
        outfile=${DIR_OUT}/v2/${statuss}/${base}
        echo -e "${filename} ${indir} ${suffix} ${outfile}" >> ${DIR_OUT}/params/v2_params.txt
done < ${tmp_v2}

echo "--- $(date) --- Create file: v2b ---"

while read file
do
        filename=$(echo -n "${file}" | cut -f 1)
        base=$(basename "$filename")
        tissue="DLPFC"

        ## Check Status ##
        if [[ "$filename" == *"Control"* ]]; then
                statuss="Control"
        else
                statuss="AD"
        fi
        indir=${DIR_FOUR}/si.results/${tissue}
        suffix="${statuss}/"
        outfile=${DIR_OUT}/v2b/${statuss}/${base}
        echo -e "${filename} ${indir} ${suffix} ${outfile}" >> ${DIR_OUT}/params/v2b_params.txt
done < ${tmp_v2b}

rm ${tmp_v2}
rm ${tmp_v2b}

echo "--- $(date) --- Done ---"
