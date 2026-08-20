#! /usr/bin/bash
#$ -cwd
#$ -N asarp
#$ -o logasarp/job.JOBNAME.JOB_ID.$TASK_ID.out
#$ -e logasarp/job.JOBNAME.JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=24:00:00,h_data=4G
#$ -pe shared 2
#$ -t 1-1:1

THREAD=2

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -e -x -o pipefail
. /u/local/Modules/default/init/modules.sh

module load picard_tools
module load samtools
module load bcftools

####################

DIR_STAR_REF=/path/to/STAR/genome/index
DIR_ASARP_CODE=/path/to/ASARP/code
DIR_PSI_CODE=/path/to/PSI_calculator/code
DIR_RDD_CODE=/path/to/RDD_pipeline/code

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/ASARP
DIR_FQ=${DIR_WORK}/FQ
DIR_BAM=${DIR_WORK}/BAM
DIR_PSI=${DIR_WORK}/PSI
DIR_RDD=${DIR_WORK}/RDD
DIR_VCF=${DIR_RDD}/VCF
DIR_CONFIG=${DIR_WORK}/CONFIG
DIR_ASARP=${DIR_WORK}/ASARP
DIR_HIST=${DIR_WORK}/HIST

ALLFQS=${DIR_FQ}/samples_names.txt
DIR_STAR=/path/to/STAR
FA=/path/to/reference/genome.fa

gene_anno=/path/to/gene/annotation/file
splice_anno=/path/to/splice/annotation/file
psi_anno=/path/to/psi/annotation/file

map=/map/individual/to/vcf
fqfilename=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${ALLFQS})
rnaseqid=$(awk -v lin=${idx} 'FNR == (lin) {print $2;}' ${ALLFQS})
indiv=$(awk -v lin=${idx} 'FNR == (lin) {print $3;}' ${ALLFQS})
vcf_sample_id=$(grep -w "${indiv}" ${map} | cut -f2 | head -n1)

echo "--- $(date) --- Download ---"

# Download BAM/fastq files from database of interest
# Conversion from BAM to fastq may be required if the database does not provide fastq files directly.

echo "--- $(date) --- Mapping ---"

${DIR_STAR}/STAR \
    --runMode alignReads \
    --runThreadN ${THREAD} \
    --genomeDir ${DIR_STAR_REF} \
    --genomeLoad NoSharedMemory \
    --readFilesIn ${DIR_FQ}/${fqfilename} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${DIR_BAM}/${rnaseqid}. \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax .04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignEndsType EndToEnd \
    --outBAMcompression 10 \
    --quantMode GeneCounts

rm -f ${DIR_FQ}/${fqfilename}

echo "--- $(date) --- Marking Duplicates ---"

java -jar ${DIR_ASARP_CODE}/scripts/MarkDuplicates.jar \
       I=${DIR_BAM}/${rnaseqid}.Aligned.sortedByCoord.out.bam \
       O=${DIR_BAM}/${rnaseqid}.rmdup.bam \
       M=${DIR_BAM}/${rnaseqid}.rmdup.metrics \
       CREATE_INDEX=true \
       VALIDATION_STRINGENCY=SILENT \
       REMOVE_DUPLICATES=true

rm -f ${DIR_BAM}/${rnaseqid}.Aligned.sortedByCoord.out.bam

echo "--- $(date) --- PSI ---"

module load python
python ${DIR_PSI_CODE}/PSI_calculator.py \
    -b ${DIR_BAM}/${rnaseqid}.rmdup.bam \
    -a ${psi_anno} \
    -p ${DIR_PSI}/${rnaseqid}.psi

echo "--- $(date) --- RDD"

chromosomes=(
     chr1
     chr2
     chr3
     chr4
     chr5
     chr6
     chr7
     chr8
     chr9
     chr10
     chr11
     chr12
     chr13
     chr14
     chr15
     chr16
     chr17
     chr18
     chr19
     chr20
     chr21
     chr22
     chrX
     chrY
)

mkdir ${DIR_RDD}/${rnaseqid}
for chrom in ${chromosomes[@]}; do
    python ${DIR_RDD_CODE}/step1_get_rdd_coordinates.py \
           -b ${DIR_BAM}/${rnaseqid}.rmdup.bam \
           --user_coordinates ${DIR_VCF}/${vcf_sample_id}.rdds.tsv \
           -c ${chrom} \
           -f ${FA} \
           -o ${DIR_RDD}/${rnaseqid}

    python ${DIR_RDD_CODE}/step2.get.mm_pileup_reads.py \
           -b ${DIR_BAM}/${rnaseqid}.rmdup.bam \
           -c ${chrom} \
           -o ${DIR_RDD}/${rnaseqid} \
           -i ${rnaseqid} \
           --bedgraph_dir ${DIR_RDD}/${rnaseqid}
done

python ${DIR_RDD_CODE}/step3.merge_reads_and_llr_cal.py \
       -o ${DIR_RDD}/${rnaseqid} \
       -i ${rnaseqid} \
       --asarp_output \
       --mono 0 \
       --mono_ratio 0.0

rm ${DIR_RDD}/${rnaseqid}.tmp*
DIR_BEDGR=${DIR_RDD}/${rnaseqid}
SNV_FILE=${DIR_RDD}/${rnaseqid}.final.all_snvs.${rnaseqid}.rdd

echo "--- $(date) --- Visualize ---"

module load R
Rscript ${DIR_ASARP_CODE}/scripts/plot.ratio_histogram.R $SNV_FILE ${rnaseqid}.allele_dist.png
mv ${rnaseqid}.allele_dist.png ${DIR_HIST}

echo "--- $(date) --- Write Config ---"

config=${DIR_CONFIG}/${rnaseqid}.config

genome_annotation=${gene_anno}
annotation_splicing_events=${splice_anno}
PSI_FILE=${DIR_PSI}/${rnaseqid}.psi

printf "snpfile\t${SNV_FILE}\n"                                > ${config}
printf "bedfolder\t${DIR_BEDGR}\n"                      >> ${config}
printf "bedext\tbed\n"                                         >> ${config}

#Annotation
printf "genefile\t${genome_annotation}\n"               >> ${config}
printf "splicingfile\t${annotation_splicing_events}\n"         >> ${config}

#Strandedness
printf "strandflag\t2\n"                                >> ${config}

#PSI params:
printf "psifile\t${PSI_FILE}\n"                                >> ${config}
printf "psitotal\t10\n"                                        >> ${config}
printf "psiexclusion\t2\n"                              >> ${config}

echo "--- $(date) --- Running ASARP ---"

perl -I ${DIR_ASARP_CODE}/asarp_vPSI ${DIR_ASARP_CODE}/asarp_vPSI/asarp.pl \
     ${DIR_ASARP}/${rnaseqid}.out_asarp \
     ${config} \
     ${DIR_ASARP_CODE}/asarp_vPSI/default.param \
     ${rnaseqid}

echo "--- $(date) --- END ---"