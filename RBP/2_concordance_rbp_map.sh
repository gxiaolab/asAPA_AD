#! /usr/bin/bash
#$ -cwd
#$ -N rosmaprbp2
#$ -o log2/job.$JOB_NAME.$TASK_ID.out
#$ -e log2/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_rt=0:15:00,h_data=1G
#$ -pe shared 1
#$ -t 1-1:1

THREAD=2

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python

####################

DIR_BASE=/home/directory/asAPA_analysis
DIR_WORK=${DIR_BASE}/RBP
DIR_ONE=${DIR_WORK}/MAPPED_step1

pval_script=./bootstrap_for_rbp_analysis.py
fdr_correction=./fdr_correction.py

DIR_OUT=${DIR_WORK}/ORGANIZE_step2
mapped_file_sig=${DIR_ONE}/asapa_snps_mapped.txt
mapped_file_ctrl=${DIR_ONE}/ctrl_snps_mapped.txt

echo "--- $(date) --- Make Content File ---"

tissue="Tissue Type" # Replace with the actual tissue name
out=${DIR_OUT}/${tissue// /_}.txt

echo "--- $(date) --- Get Tissue and Status ---"

tmp=${DIR_OUT}/${tissue// /_}_tmp.txt
tmp_N=${DIR_OUT}/${tissue// /_}_nonsig_tmp.txt

grep "${tissue}" ${mapped_file_sig} >> ${tmp}
grep "${tissue}" ${mapped_file_ctrl} >> ${tmp_N}

echo "--- $(date) --- Get Tuples  ---"

tmp_2=${DIR_OUT}/${tissue// /_}_tmp2.txt
tmp_2_N=${DIR_OUT}/${tissue// /_}_nonsig_tmp2.txt
tmp_3=${DIR_OUT}/${tissue// /_}_tmp3.txt
tmp_3_N=${DIR_OUT}/${tissue// /_}_nonsig_tmp3.txt

awk -F'\t' -v OFS='\t' '{print $2, $8}' ${tmp} >> ${tmp_2}
awk -F'\t' -v OFS='\t' '{print $2, $8}' ${tmp_N} >> ${tmp_2_N}

sort -u ${tmp_2} >> ${tmp_3}
sort -u ${tmp_2_N} >> ${tmp_3_N}

echo "--- $(date) --- Count  ---"

tmp_4=${DIR_OUT}/${tissue// /_}_tmp4.txt

tiss_gene_count=$(awk -F'\t' '{print $1}' ${tmp_3} | sort -u | wc -l)
tiss_gene_count_N=$(awk -F'\t' '{print $1}' ${tmp_3_N} | sort -u | wc -l)

while read pair
do
        rbp=$(echo -n "${pair}" | cut -f 2 | tr -d $'\n')
        test_count=$(cut -f2 ${tmp_3} | grep -w "${rbp}" | wc -l || true)
        ctrl_count=$(cut -f2 ${tmp_3_N} | grep -w "${rbp}" | wc -l || true)
        echo -e "${rbp}\t${test_count}\t${ctrl_count}\t${tiss_gene_count}\t${tiss_gene_count_N}" >> ${tmp_4}
done < ${tmp_3}

sort -u ${tmp_4} >> ${out}

rm ${tmp}
rm ${tmp_N}
rm ${tmp_2}
rm ${tmp_2_N}
rm ${tmp_3}
rm ${tmp_3_N}
rm ${tmp_4}

echo "--- $(date) --- Done  ---"
