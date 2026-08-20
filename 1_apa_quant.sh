#! /usr/bin/bash
#$ -cwd
#$ -N apa_quant
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=24:00:00,h_data=16G
#$ -pe shared 4
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python
module load bedtools

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_OUT=${DIR_WORK}/aQTL

gtf=/path/to/annotation/gtf/file
polyadb=/path/to/polyA_DB2/reference/reference.bed  # File with simplified coordinates of proximal and distal APA regions of interest for each gene (chrom | proximal lower bound | proximal upper bound | distal lower bound | distal upper bound | strand)
polyadb_ref=/path/to/polyA_DB/human.PAS.map_withPAS_on_Refseq_NM.txt    # File with full information from polyADB database

####################

echo "--- $(date) --- Get events from ${gtf} and sort ---"

mkdir -p ${DIR_OUT}

python ./apa_quant1.py \
-i ${gtf} \
-o ${DIR_OUT}/asapa_events.txt \
-b ${DIR_OUT}/asapa_sites.bed

awk '{
    split($4,a,"_"); 
    gene=a[1]; 
    count[gene]++; 
    lines[NR]=$0; 
    genes[NR]=gene
}
END{
    for(i=1;i<=NR;i++){
        if(count[genes[i]]>1){
            print lines[i]
        }
    }
}' ${DIR_OUT}/asapa_sites.bed > ${DIR_OUT}/asapa_sites_filtered.bed

echo "--- $(date) --- Get events from ${polyadb} and sort ---"

python ./apa_quant2.py ${polyadb} ${polyadb_ref} ${DIR_OUT}/polyadb_sites.bed

awk '{
    split($4,a,"_"); 
    gene=a[1]; 
    count[gene]++; 
    lines[NR]=$0; 
    genes[NR]=gene
}
END{
    for(i=1;i<=NR;i++){
        if(count[genes[i]]>1){
            print lines[i]
        }
    }
}' ${DIR_OUT}/polyadb_sites.bed > ${DIR_OUT}/polyadb_sites_filtered.bed

echo "--- $(date) --- Get coverage of bed files ---"

bams=${DIR_OUT}/bam_list.txt    # List of all bam names

bedtools multicov -bams $(paste -sd' ' ${bams}) -bed ${DIR_OUT}/asapa_sites_filtered.bed > ${DIR_OUT}/bed_coverage/asapa_raw_counts.txt
bedtools multicov -bams $(paste -sd' ' ${bams}) -bed ${DIR_OUT}/polyadb_sites_filtered.bed > ${DIR_OUT}/bed_coverage/polyadb_raw_counts.txt

python ./apa_quant3.py \
    ${DIR_OUT}/bed_coverage/asapa_raw_counts.txt \
    ${DIR_OUT}/bam_list.txt \
    ${DIR_OUT}/matrices/asapa_raw_matrix.txt \
    ${DIR_OUT}/matrices/asapa_relative_usage.txt

python ./apa_quant3.py \
    ${DIR_OUT}/bed_coverage/polyadb_raw_counts.txt \
    ${DIR_OUT}/bam_list.txt \
    ${DIR_OUT}/matrices/polyadb_raw_matrix.txt \
    ${DIR_OUT}/matrices/polyadb_relative_usage.txt

echo "--- $(date) --- Done ---"
