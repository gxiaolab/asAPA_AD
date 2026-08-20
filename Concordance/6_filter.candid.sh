#! /usr/bin/bash
#$ -cwd
#$ -N filter.candid
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

echo "--- $(date) --- Set vars ---"

# if rerunning, just do bedtools intersect, uniq with 3utr file instead of this script

candid_files=(/path/to/Concordance/STEP3_VCFclean/candid.pass.common.nonegative.bed/*)
infile=${candid_files[$idx-1]}
id=$(basename $infile | cut -d. -f1)
OUT_DIR=/path/to/Concordance/STEP3_VCFclean/candid.bed.filter
outfile=${OUT_DIR}/${id}.candid.filter.bed
utr_bed=/path/to/3UTR/annotation.bed
echo Input:$infile
echo Output:$outfile

echo "--- $(date) --- Intersect ---"

bedtools intersect -a $utr_bed -b $infile -wb \
        | awk -v OFS='\t' '{print $13,$14,$15,$16}' > $outfile

echo "--- $(date) --- Done ---"
