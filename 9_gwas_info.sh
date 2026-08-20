#! /usr/bin/bash
#$ -cwd
#$ -N cleangwas
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=1:00:00,highp
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_GWAS=${DIR_WORK}/GWAS
DIR_INFO=${DIR_GWAS}/GWAS_SNPs
DIR_STEP0=${DIR_GWAS}/STEP0_CLEANINFO

gwas_files=${DIR_INFO}/gwas_filenames.txt

echo "--- $(date) --- Clean ---"

while read gwas_file
do
    gwas=${DIR_INFO}/${gwas_file}
    disease="${gwas_file%.tsv}"
    out=${DIR_STEP0}/${disease}.txt
    tmp1=${DIR_STEP0}/${disease}_col1_tmp.txt
    tmp13=${DIR_STEP0}/${disease}_col13_tmp.txt
    tmp_both=${DIR_STEP0}/${disease}_1_13_tmp.txt
    cut -f1 ${gwas} | grep "chr" | grep -v "locations" | sed -E 's/(-).*$//' | tr -d '"' | tr ',' '\n' | tr ':' '\t' >> ${tmp1} || true
    cut -f13 ${gwas} | awk -F'\t' '$1!="-"' | grep -v "locations" | tr -d '"' | tr ',' '\n' | sed 's/^/chr/' | tr ':' '\t' >> ${tmp13} || true
    cat ${tmp1} ${tmp13} >> ${tmp_both}
    sort -u ${tmp_both} >> ${out}

    rm ${tmp1}
    rm ${tmp13}
    rm ${tmp_both}
done < ${gwas_files}

echo "--- $(date) --- Done ---"
