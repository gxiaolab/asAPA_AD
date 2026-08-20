#! /usr/bin/bash
#$ -cwd
#$ -N si_scores
#$ -o log2/job.$JOB_NAME.$TASK_ID.out
#$ -e log2/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_data=2G,h_rt=2:00:00
#$ -pe shared 2
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
DIR_ONE=${DIR_CONC}/STEP1_PREPROC
DIR_TWO=${DIR_CONC}/STEP2_RDD
DIR_OUT=${DIR_CONC}/STEP4_SISCORES

genotype_metadata=/path/to/genotype/metadata.txt
metadata=${DIR_WORK}/ASARP/FQ/samples_names.txt

###################

echo "--- $(date) --- Set variables ---"

tag_files=(${DIR_TWO}/*.tag.rdd.bed)
tag=${tag_files[$idx-1]}
rnaseqid=$(basename $tag | cut -d "." -f 1)

indv_id=$(awk -v rnaseq="$rnaseqid" -F "\t" '$2==(rnaseq) {print $3}' ${metadata})
echo individual id: $indv_id

vcf_sample_id=$(less ${genotype_metadata} | grep $indv_id | awk '{print $2}')
echo vcf sample id: $vcf_sample_id

# Previously said .candid.common.bed. Why?
candid=/path/to/analysis/Concordance/STEP3_VCFclean/candid.bed.filter/${vcf_sample_id}.candid.filter.bed

annot=/path/to/analysis/Concordance/gencode.v38.GRCh38.bed
maxD=0.5
search=INF

het=/path/to/analysis/Concordance/STEP3_VCFclean/tag.bed.heterozygous/${vcf_sample_id}.het.bed

outfile=${DIR_OUT}/${rnaseqid}.splicing.out

echo "--- $(date) --- Calc SI ---"

python ./splicing.concordance.tag.het.py -i $candid -t $tag -d $rnaseqid -m $maxD -s $search -a $annot -o $outfile -b $het

echo "--- $(date) --- Organize ---"

tissue="Tissue Type" #replace with tissue type
out_ctrlad=${DIR_OUT}/si.results
name=$(basename $outfile)
mkdir -p ${DIR_OUT}/by_chrom

mapping=${DIR_WORK}/path/to/AD_status/metadata.txt

for i in {1..22} #X Y #no X or Y!!!
do
    less $outfile | grep -w chr$i > ${DIR_OUT}/by_chrom/chr${i}.$name || true
    status=$(less $mapping | grep $indv_id | awk '{print $2}')
    mkdir -p $out_ctrlad/$tissue/$status/
    ln -s ${DIR_OUT}/by_chrom/chr${i}.$name $out_ctrlad/$tissue/$status/
done

echo "--- $(date) --- Done ---"
