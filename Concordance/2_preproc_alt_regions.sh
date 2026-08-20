#! /usr/bin/bash
#$ -cwd
#$ -N preproc_reg
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
DIR_OUT=${DIR_CONC}/STEP1_PREPROC

snpfull=${DIR_OUT}/universal_tag_snp_uniq.bed
annotscript=/path/to/ASARP_code/annotSnvGene.pl
genepred=/path/to/gtf/file.gtf
config=/path/to/ASARP/config/files/sample1.config
extractscript=/path/to/concordance/scripts/universal.extractregions.py

####################

echo "--- $(date) --- Make Main ---"

module load perl

rm -f ${DIR_OUT}/tmp.bed
rm -f ${DIR_OUT}/extractregions.log
rm -f ${DIR_OUT}/universal.tag.regions.bed
rm -f ${DIR_OUT}/universal.tag.regionsfinal.bed

snps=${DIR_OUT}/tmp.bed
snpfullextra=${DIR_OUT}/universal_tag_snp_uniq_extra.bed
awk -F'\t' -v OFS=" " '{print $1, $2 + 1, $3}' ${snpfullextra} > ${snps}

perl -I /path/to/ASARP_code/ASARP ${annotscript} $snps $genepred $config > ${DIR_OUT}/extractregions.log

while IFS=$'\t' read -r chr pos strd gene x type ref alt; do
    posaddone=$((pos + 1))

    region=$(awk -F',' -v p="$posaddone" '$2 == p {print $3; exit}' ${DIR_OUT}/extractregions.log)

    if [[ -z "$region" || "$region" != *-* ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$chr" "$pos" "$strd" "$gene" "$region" "$type" "$ref" "$alt" \
            >> ${DIR_OUT}/universal.tag.regions.missing.bed
        continue
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$chr" "$pos" "$strd" "$gene" "$region" "$type" "$ref" "$alt" \
        >> ${DIR_OUT}/universal.tag.regions.bed
done < ${snpfull}

python3 ${extractscript} ${DIR_OUT}/universal.tag.regions.bed ${DIR_OUT}/universal.tag.regionsfinal.bed

less ${DIR_OUT}/universal.tag.regionsfinal.bed | cut -d';' -f1 > ${DIR_OUT}/universal.tag.regionsfinal_clean.bed

echo "--- $(date) --- Done ---"
