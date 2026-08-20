#! /usr/bin/bash
#$ -cwd
#$ -N rddtag
#$ -o log/job.$JOB_NAME.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_rt=12:00:00,h_data=16G
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh

module load samtools
module load anaconda3
module load python

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_CONC=${DIR_WORK}/Concordance
DIR_ASARP=${DIR_WORK}/ASARP/asarp_final
DIR_ASARPOUT=${DIR_WORK}/ASARP/ASARP
events=${DIR_ASARP}/path/to/significant/asAPA_events.txt
DIR_ONE=${DIR_CONC}/STEP1_PREPROC
DIR_OUT=${DIR_CONC}/STEP2_RDD

snpfull=${DIR_ONE}/universal_tag_snp_rdd_clean.bed
rnaseqids=${DIR_WORK}/ASARP/FQ/samples_names.txt

DIR_RDD_PIPELN=/path/to/RDD_scripts
genome_fa=/path/to/genome/hg38.fa

####################

echo "--- $(date) --- Set variables ---"

DIR_BAM=${DIR_WORK}/ASARP/BAM
rnaseqid=$(awk -v lin=${idx} 'FNR == (lin) {print $2;}' ${rnaseqids})
bamfile=${DIR_BAM}/${rnaseqid}.rmdup.bam

echo "--- $(date) --- RDD ---"

chroms=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
        chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)

for chrom in "${chroms[@]}"
do
        python ${DIR_RDD_PIPELN}/step1_get_rdd_coordinates.py \
                -b $bamfile \
                --user_coordinates ${snpfull} \
                -c $chrom \
                -o ${DIR_OUT}/${rnaseqid} \
                -f ${genome_fa} \

        python ${DIR_RDD_PIPELN}/step2.get.mm_pileup_reads.py \
               -b $bamfile \
               -c $chrom \
               -o ${DIR_OUT}/${rnaseqid} \
               -i ${rnaseqid} 
done

echo $(date) step 3: merge reads, calculate LLR
python ${DIR_RDD_PIPELN}/step3.merge_reads_and_llr_cal.py \
       -o ${DIR_OUT}/${rnaseqid} \
       -i $rnaseqid \
       --asarp_output \
       --mono 0 \
       --mono_ratio 0.0

echo $(date) remove temporary files
rm ${DIR_OUT}/${rnaseqid}.tmp*

echo "--- $(date) --- Done ---"
