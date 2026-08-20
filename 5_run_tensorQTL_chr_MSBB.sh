#!/bin/bash
#$ -cwd
#$ -o log/tensor_chr_MSBB.$JOB_ID.$TASK_ID
#$ -j y
#$ -l h_rt=12:00:00,h_data=36G
#$ -pe shared 1
#$ -t 1-88:1

source ~/.bash_profile 
source ~/miniconda3/bin/activate plink   

##### 1. cis mode 
pheno=asapa_relative_usage.rename
chr=`cat run.tensorQTL.chr.MSBB | awk '{print $1}' | sed -n ${SGE_TASK_ID}p`
region=`cat run.tensorQTL.chr.MSBB | awk '{print $2}' | sed -n ${SGE_TASK_ID}p`

grep $chr -w ./${pheno}.${region}.normalize.wgs_samples.bed | sort -k1,1V -k2,2n -k3,3n | cat <(head -n1 ./${pheno}.${region}.normalize.wgs_samples.bed) - > 02_phenotype/chr/${pheno}.${region}.normalize.wgs_samples.sort.${chr}.bed 

plink_prefix_path=./all_chr.final_rename_${region}.${chr}
pheno_bed=./chr/${pheno}.${region}.normalize.wgs_samples.sort.${chr}.bed
covariates_file=./${pheno}.${region}.normalize.wgs_samples.braak.covariants
out_prefix=03_tensorQTL_braak_out/cis_res/${pheno}.${region}.${chr} 

test -d 03_tensorQTL_braak_out/cis_res/ || mkdir -p 03_tensorQTL_braak_out/cis_res/ 

### chromosome: chr1, chr2 ... 
### the order of samples in pheno_bed and covariants_file should be same 
### the value of covariantes should be number  
python3 -m tensorqtl ${plink_prefix_path} ${pheno_bed} ${out_prefix} \
    --covariates ${covariates_file} \
    --mode cis  \
    --window 1000000 && echo finish_QTL $pheno $chr $region



