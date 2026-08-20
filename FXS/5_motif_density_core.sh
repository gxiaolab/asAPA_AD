#! /usr/bin/bash
#$ -cwd
#$ -N mot_dens_core
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=2:00:00,h_data=2G,highp
#$ -pe shared 1
#$ -t 1-5:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_FX=${DIR_WORK}/FXS
DIR_FMRP=${DIR_FX}/FMRP_eCLIP
DIR_MOTIF=${DIR_FX}/FMRP_Motif
DIR_PDB=${DIR_FMRP}/polyA_DB

DIR_FOUR=${DIR_PDB}/4_coverage
DIR_SIX=${DIR_PDB}/6_usage
info=${DIR_FX}/info/sample_to_disease_status.txt
cohorts=${DIR_FOUR}/cohorts.txt

genome_fasta=/path/to/genome/fasta.fa
gene_map=/path/to/map/geneid_to_genename.txt
gene_pred=/path/to/gtf_file
map=${DIR_PDB}/human.PAS.map_withPAS_on_Refseq_NM.txt

echo "--- $(date) --- Clean ---"

cohort=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${cohorts})
region_file=${DIR_MOTIF}/density_core/${cohort}_regions_case_utr_txsearch.txt

tmp=${DIR_MOTIF}/density_core/${cohort}_tmp.txt
clean_file=${DIR_MOTIF}/density_core/${cohort}_regions_case_utr_txsearch_clean.txt

awk 'NF==10 && !/^[ \t]*$/ && !($0 ~ /[[:space:]]{2,}/)' ${region_file} | awk -F'\t' '$5 - $4 >= 400' | awk -F'\t' '$5 - $4 < 20000' | awk -F'\t' '$8==$9' | awk -F'\t' '!seen[$4,$5,$7,$8]++' >> ${tmp}
(echo -e "chrom\tgene\tPAS\tregion_start\tregion_end\ttxend\tevent_category\tstrand\tstrand_gencode\tcase"; cat ${tmp}) >> ${clean_file}

rm ${tmp}

echo "--- $(date) --- Scan Motifs ---"

window=500
slide=50
smooths=(50)
for smooth in "${smooths[@]}"
do
  out=${DIR_MOTIF}/density_core/${cohort}_${window}_${smooth}_final_ext
  /u/home/k/kofiamoa/.conda/envs/ML/bin/python ./motif_density_core.py \
    -i ${clean_file} \
    -g ${genome_fasta} \
    -o ${out} \
    -w ${window} \
    --sliding ${slide} \
    --smooth ${smooth}
done

rm ${clean_file}

echo "--- $(date) --- Done ---"
