#! /usr/bin/bash
#$ -cwd
#$ -N polyadb_diff
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=1:00:00,h_data=4G,highp
#$ -pe shared 1
#$ -t 1-1:1

THREAD=2

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python
module load bedtools

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_FX=${DIR_WORK}/FXS
DIR_FMRP=${DIR_FX}/FMRP_eCLIP
DIR_WIG=${DIR_FX}/WIG
DIR_PDB=${DIR_FMRP}/polyA_DB
pas_all=${DIR_OUT}/pas.bed
endt_all=${DIR_OUT}/transcript_ends.bed
names=${DIR_WIG}/wig_names.txt

DIR_FOUR=${DIR_OUT}/4_coverage
info=${DIR_FX}/info/map_sample_to_fxs_status.txt
cohorts=${DIR_FOUR}/cohorts.txt

DIR_FIVE=${DIR_OUT}/5_average
average=./polyA_DB_diffAPA_average.py

DIR_SIX=${DIR_OUT}/6_usage
usage=./polyA_DB_diffAPA_usage.py

combine=./polyA_DB_diffAPA_combine.py

echo "--- $(date) --- Sort polyDB info ---"

bedgraph=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${names})
sraid=${bedgraph%.wig}

bedtools sort -i ${pas_all} > ${DIR_OUT}/pas_sorted.bed
bedtools sort -i ${endt_all} > ${DIR_OUT}/transcript_ends_sorted.bed

echo "--- $(date) --- Get Wig, Sort, and Map ---"

bedtools sort -i ${DIR_WIG}/${bedgraph} > ${DIR_FOUR}/${sraid}_sorted.bed

cohort=$(grep ${sraid} ${info} | cut -f4)
stat=$(grep ${sraid} ${info} | cut -f2)

mkdir -p ${DIR_FOUR}/${cohort}

bedtools map -a ${DIR_OUT}/pas_sorted.bed -b ${DIR_FOUR}/${sraid}_sorted.bed -c 4 -o mean > ${DIR_FOUR}/${cohort}/${sraid}_${stat}_PAS_coverage.bed
bedtools map -a ${DIR_OUT}/transcript_ends_sorted.bed -b ${DIR_FOUR}/${sraid}_sorted.bed -c 4 -o mean > ${DIR_FOUR}/${cohort}/${sraid}_${stat}_transcript_coverage.bed

rm ${DIR_FOUR}/${sraid}_sorted.bed

echo "--- $(date) --- Average the coverage for each cohort & Status ---"

cohort=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${cohorts})

cohort_to_average=${DIR_FOUR}/${cohort}
out_path=${DIR_FIVE}/${cohort}

mkdir -p ${out_path}

python ${average} ${cohort_to_average} ${out_path}

echo "--- $(date) --- Get Relative Usage metric ---"

cohort=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${cohorts})

# Reference is built from all proximal and distal apa regions
ref=${DIR_OUT}/reference/reference.bed
fxs_pas=${DIR_FIVE}/${cohort}/FXS_avg_PAS_coverage.bed
ctrl_pas=${DIR_FIVE}/${cohort}/Control_avg_PAS_coverage.bed
fxs_transcript=${DIR_FIVE}/${cohort}/FXS_avg_transcript_coverage.bed
ctrl_transcript=${DIR_FIVE}/${cohort}/Control_avg_transcript_coverage.bed
out_usage=${DIR_SIX}/${cohort}/all_usage_difference_.5.txt

mkdir -p ${DIR_SIX}/${cohort}

python ${usage} ${ref} ${fxs_pas} ${ctrl_pas} ${fxs_transcript} ${ctrl_transcript} ${out_usage}

final_usage=${DIR_SIX}/${cohort}/final_usage_difference_.5.txt
grep -v 'NA' ${out_usage} > ${final_usage}

echo "--- $(date) --- Done ---"
