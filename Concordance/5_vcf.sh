#! /usr/bin/bash
#$ -cwd
#$ -N isec_vcf
#$ -o log3/job.$JOB_NAME.$TASK_ID.out
#$ -e log3/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_data=2G,h_rt=24:00:00 # 2G
#$ -pe shared 4 # 4
#$ -t 1-1:1

idx=$SGE_TASK_ID
THREAD=4

source $HOME/.bash_profile

set -e -x -o pipefail
. /u/local/Modules/default/init/modules.sh
module load anaconda3
module load bedtools

###################

DIR_BASE=/home/directory
DIR_THREE=${DIR_BASE}/Concordance/STEP3_VCFclean
bcftools=/path/to/bcftools

vcfs=${DIR_BASE}/path/to/VCFs/sample_names.txt

###################

echo "--- $(date) --- Prep Merged VCF ---"

# RUN THIS JUST ONCE WITH TASK ID 1-1
merged_vcf=${DIR_BASE}/VCF/hg38/merged/genotype_merged.vcf.gz
${bcftools} view -Oz -o ${DIR_BASE}/VCF/hg38/merged/genotype_merged.bgzip.vcf.gz ${merged_vcf}
${bcftools} index --tbi --threads 4 ${DIR_BASE}/VCF/hg38/merged/genotype_merged.bgzip.vcf.gz

merged_vcf=${DIR_BASE}/VCF/hg38/merged/genotype_merged.bgzip.vcf.gz

echo "--- $(date) --- Filter VCFs ---"

id=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${vcfs})
DIR_FILT=${DIR_BASE}/path/to/VCFs/filtered_vcfs

sample_vcf=${DIR_FILT}/${id}.filtered.vcf.gz

if [[ ! -f ${sample_vcf} ]]
then
${bcftools} view \
    --threads ${THREAD} \
    --samples ${id} \
    --force-samples \
    --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
    --output-type z \
    ${merged_vcf} | \
${bcftools} filter \
    --threads ${THREAD} \
    -i 'FILTER="PASS"' \
    --output-type z | \
${bcftools} view \
    --threads ${THREAD} \
    --types snps \
    -m2 -M2 \
    --output-type z | \
${bcftools} annotate \
    --threads ${THREAD} \
    --set-id '%CHROM:%POS' \
    --output-type z \
    --output ${sample_vcf}

${bcftools} index --tbi --threads ${THREAD} ${sample_vcf}
fi

echo "--- $(date) --- dbSNP ---"

infile=${sample_vcf}
OUT_DIR=${DIR_THREE}/candid.pass.common.nonegative.bed
common_dbsnp=/path/to/dbSNP/dbsnp.38.hg38.common.vcf.gz
if [[ -f ${OUT_DIR}/${id}.candid.common.bed ]]
then
    echo "This ID is already complete."
    exit 0
fi
bedtools intersect -header -u -a ${infile} -b ${common_dbsnp} > ${OUT_DIR}/${id}.filtered.common.vcf

echo "--- $(date) --- Get candid bed ---"

# output in candid.bed format
${bcftools} query -i 'TYPE="snp" && GT!="./." && NEGATIVE_TRAIN_SITE!=1' -f'%CHROM\t%POS0\t%END\t[%GT]|GT|%REF>%ALT{0}\n' \
                    --output ${OUT_DIR}/${id}.candid.common.bed \
                    ${OUT_DIR}/${id}.filtered.common.vcf

rm ${OUT_DIR}/${id}.filtered.common.vcf

echo "--- $(date) --- Done ---"
