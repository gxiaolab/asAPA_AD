#! /usr/bin/bash
#$ -cwd
#$ -N preproc_tag
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=4:00:00,h_data=4G
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_CONC=${DIR_WORK}/Concordance
DIR_ASARP=${DIR_WORK}/ASARP/asarp_final
DIR_ASARPOUT=${DIR_WORK}/ASARP/ASARP
events=${DIR_ASARP}/path/to/significant/asAPA_events.txt
DIR_OUT=${DIR_CONC}/STEP1_PREPROC

####################

echo "--- $(date) --- Make Main ---"

rm -f ${out} ${out2} ${out3}

out=${DIR_OUT}/universal_tag_snp_all.bed

while read event
do
        rnaseqid=$(echo -n "${event}" | cut -f 1)
        gene=$(echo -n "${event}" | cut -f 3)
        chrom=$(echo -n "${event}" | cut -f 4)
        pos=$(echo -n "${event}" | cut -f 5)
        pos=$((pos - 1))
        snpinfo=pos=$(echo -n "${event}" | cut -f 13)
        reference=$(echo -n "${snpinfo}" | cut -d',' -f3 | cut -d'>' -f1)
        alternate=$(echo -n "${snpinfo}" | cut -d',' -f3 | cut -d'>' -f2)
        strand=$(echo -n "${event}" | cut -f 24)

        echo -e "${chrom}\t${pos}\t${strand}\t${gene}\t0:0\taltterm\t${reference}\t${alternate}" >> ${out}
done < ${events}

echo "--- $(date) --- Get uniq ---"

out2=${DIR_OUT}/universal_tag_snp_uniq.bed

sort -u ${out} >> ${out2}

echo "--- $(date) --- Get RDD ---"

out3=${DIR_OUT}/universal_tag_snp_rdd.bed

cut -f1,2,7,8 ${out2} >> ${out3}

echo "--- $(date) --- Done ---"
