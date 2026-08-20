#! /usr/bin/bash
#$ -cwd
#$ -N deepripeB
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=1:00:00,highp
#$ -pe shared 1
#$ -t 1-1:1

THREAD=2

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_OUT=${DIR_WORK}/RBP/DeepRiPe
DIR_TOOLS=/path/to/DeepRiPe/scripts
DIR_CONC=${DIR_WORK}/Concordance/CONCORDANCE_RERUN_CORRECTED
DIR_VAR=${DIR_OUT}/variants
conc_events=${DIR_CONC}/path/to/putative_functional_snps.txt
out=${DIR_VAR}/snps_corrected.bed

####################

echo "--- $(date) --- Run ---"

while read event
do
        chrom=$(echo -n "${event}" | cut -f 1)
        pos=$(echo -n "${event}" | cut -f 2)
        gene=$(echo -n "${event}" | cut -f 20)
        #nev=$(echo -n "${event}" | cut -f 7)
        strand=$(echo -n "${event}" | cut -f 6)
        tissue=$(echo -n "${event}" | cut -f 7)
        tissue_1="${tissue#[}"
        tissue_2="${tissue_1%]}"
        tissue_clean="$(echo "${tissue_2}" | tr ' ' '_')"
        statuss=$(echo -n "${event}" | cut -f 8)
        reference=$(echo -n "${event}" | cut -f 17)
        alternate=$(echo -n "${event}" | cut -f 18)
        annot=$(echo -n "${event}" | cut -f 19)
        sample_size=$(echo -n "${event}" | cut -f 9)

        if (( ${sample_size} >= 100 )); then
                confidence="High_confidence"
        elif (( ${sample_size} < 50 )); then
                confidence="Low_confidence"
        else
                confidence="Medium_confidence"
        fi

        echo -e "${chrom}\t${pos}\t${pos}\t${annot}|${statuss}|${chrom}:${pos}|${gene}|${strand}|alternate|${alternate}|${confidence}\t${alternate}/${reference}\t${strand}" >> ${out}

done < ${conc_events}

echo "--- $(date) --- Done ---"
