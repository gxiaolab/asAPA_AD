#! /usr/bin/bash
#$ -cwd
#$ -N polyadb_type
#$ -o log2/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log2/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=24:00:00,h_data=4G
#$ -pe shared 2
#$ -t 1-22:1

THREAD=2

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_FX=${DIR_WORK}/FXS
DIR_FMRP=${DIR_FX}/FMRP_eCLIP
DIR_PDB=${DIR_FMRP}/polyA_DB

mapped=${DIR_PDB}/human.PAS.map_withPAS_on_Refseq_NM.txt
polyadb_data=${DIR_PDB}/human.PAS.txt

echo "--- $(date) --- Get info ---"

out=${DIR_PDB}/by_chrom/chr${idx}_type_annotated_human.PAS_on_Refseq.txt
tmp=${DIR_PDB}/by_chrom/chr${idx}.tmp.txt

grep -w "chr${idx}" ${mapped} >> ${tmp}

while read event
do
        line=$(echo -n "${event}" | tr -d ' ')
        gene=$(echo -n "${event}" | cut -f6 | tr -d ' ')
        pasid=$(echo -n "${event}" | cut -f7 | tr -d ' ')
        info_length=$(grep -w "${pasid}" ${polyadb_data} | awk -F'\t' -v gn="${gene}" '$9==gn' | wc -l || true)
        if [ ${info_length} == "1" ]; then
                pastype=$(grep -w "${pasid}" ${polyadb_data} | awk -F'\t' -v gn="${gene}" '$9==gn' | cut -f14 | tr -d ' ')
                echo -e "${line}\t${pastype}" >> ${out}
        else
                echo -e "More than one entry with this PAS ID: ${pasid} - ${info_length}!"
        fi
done < ${tmp}

rm ${tmp}

echo "--- $(date) --- Get 3'UTR APA ---"

out2=${DIR_PDB}/2_3utr_by_chrom/chr${idx}_3utr_annotated_human.PAS_on_Refseq.txt

grep "3'UTR" ${out} | awk -F'\t' '!seen[$2,$10]++' >> ${out2}

echo "--- $(date) --- Done ---"
