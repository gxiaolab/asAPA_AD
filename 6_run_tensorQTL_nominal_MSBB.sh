#!/bin/bash
#$ -cwd
#$ -o log/nominal_MSBB.$JOB_ID.$TASK_ID
#$ -j y
#$ -l h_rt=6:00:00,h_data=36G
#$ -pe shared 1
#$ -t 1-4:1

source ~/.bash_profile 
source ~/miniconda3/bin/activate plink 

#### output all tested SNPs for downstream coloc analysis 
test -d 03_tensorQTL_braak_out/nominal_res/ || mkdir -p 03_tensorQTL_braak_out/nominal_res/ 
pheno=asapa_relative_usage.rename
region=`cat run.tensorQTL.chr.MSBB | awk '{print $2}' | sort -u | sed -n ${SGE_TASK_ID}p`

sort -k1,1V -k2,2n -k3,3n <(tail -n +2 ./${pheno}.${region}.normalize.wgs_samples.bed) | cat <(head -n1 ./${pheno}.${region}.normalize.wgs_samples.bed) - > 02_phenotype/${pheno}.${region}.normalize.wgs_samples.sort.bed

plink_prefix_path=./all_chr.final_rename_${region}
pheno_bed=./${pheno}.${region}.normalize.wgs_samples.sort.bed
covariates_file=./${pheno}.${region}.normalize.wgs_samples.braak.covariants
out_prefix=03_tensorQTL_braak_out/nominal_res/${pheno}.${region}.nominal

python3 -m tensorqtl ${plink_prefix_path} ${pheno_bed} ${out_prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal \
    --window 1000000


