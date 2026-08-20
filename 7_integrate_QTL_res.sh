#!/bin/bash
#$ -cwd
#$ -o log/integrate.$JOB_ID.$TASK_ID
#$ -j y
#$ -l h_rt=2:00:00,h_data=36G
#$ -pe shared 1
#$ -t 1-4:1

source ~/.bash_profile 
source ~/miniconda3/bin/activate plink 

test -d 03_tensorQTL_braak_out/cis_res_integrate/ || mkdir 03_tensorQTL_braak_out/cis_res_integrate/ 

pheno=asapa_relative_usage.rename
region=`sed -n ${SGE_TASK_ID}p <(awk '{print $2}' run.tensorQTL.chr.MSBB | sort -u)` 

cis_out_prefix="03_tensorQTL_braak_out/cis_res/$pheno.$region."
nominal_out_prefix="03_tensorQTL_braak_out/nominal_res/$pheno.$region.nominal"
integrate_cis_out="03_tensorQTL_braak_out/cis_res_integrate/$pheno.$region.top.cis_qtl.txt"
allpair_fdr_out="03_tensorQTL_braak_out/cis_res_integrate/$pheno.$region.all.cis_qtl.txt"

python3 f2_integrate_QTL.py $cis_out_prefix $nominal_out_prefix $integrate_cis_out $allpair_fdr_out

echo "output: $integrate_cis_out"
echo "output: $allpair_fdr_out"




