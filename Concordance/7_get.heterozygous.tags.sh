#! /usr/bin/bash
#$ -cwd
#$ -N get.het.tags
#$ -o log2/job.$JOB_NAME.$TASK_ID.out
#$ -e log2/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_data=8G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -e -x -o pipefail
. /u/local/Modules/default/init/modules.sh
module load bedtools

###################

echo "--- $(date) --- Intersect ---"

DIR_OUT=/path/to/Concordance/STEP3_VCFclean/tag.bed.heterozygous
univ_tag=/path/to/Concordance/STEP1_PREPROC/tag.bed/universal.tag.uniq.overlapformat.bed

for candid in /path/to/Concordance/STEP3_VCFclean/candid.bed.filter/*; do
    id=$(basename "$candid" | cut -d "." -f 1)

    echo "Processing ${id}"

    bedtools intersect -a "$candid" -b "$univ_tag" | grep "0/1" > "${DIR_OUT}/${id}.het.bed" || true
done

echo "--- $(date) --- Done ---"
