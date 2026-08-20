#!/bin/bash
#$ -cwd
#$ -o log/prepare_MSBB.$JOB_ID
#$ -j y
#$ -l h_rt=03:00:00,h_data=20G
#$ -pe shared 1

source ~/.bash_profile 
source ~/miniconda3/bin/activate qtl  

for region in FP IFG PHG STG 
do 
    meta_path=./specimenID_syn_region_WGS_raceW.txt
    genotype_PC_path=./all_chr.final_rename_${region}.PCs.eigenvec 
    stats_path=./MSBB_${region}/seq_stat.mat
    mat_path=./asapa_relative_usage.rename.${region}.filter.txt

    Rscript f1_prepare_QTLinput_MSBB.R $mat_path $meta_path $genotype_PC_path asapa_relative_usage.rename.${region} $stats_path $region && echo finish $region   
done 








