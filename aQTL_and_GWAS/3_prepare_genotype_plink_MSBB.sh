#!/bin/bash
#$ -cwd
#$ -o log/prepare_plink_MSBB.$JOB_ID
#$ -j y
#$ -l h_rt=2:00:00,h_data=10G 
#$ -pe shared 1 

source ~/.bash_profile 
source ~/miniconda3/bin/activate plink 

meta=./specimenID_syn_region_WGS_raceW.txt  # metadata with rnaseqid to brain region

for region in FP IFG PHG STG 
do 
	grep $region $meta | awk 'OFS="\t"{print $5,$5,$1,$1}' > ./${region}_rename_id

    plink --bfile ./all_chr.final --keep <(cut -f 1,2 ./${region}_rename_id) --make-bed --keep-allele-order --out ./all_chr.final.${region} 
    plink --bfile ./all_chr.final.${region} --update-ids ./${region}_rename_id --make-bed --keep-allele-order --out ./all_chr.final_rename_${region}
    sed -i 's/^/chr/' ./all_chr.final_rename_${region}.bim

    rm ./all_chr.final.${region}.* 

    ##### separate chromosome 
    for i in {1..22}
    do
        chr="chr${i}"
        plink --bfile ./all_chr.final_rename_${region} --chr $chr --allow-extra-chr --make-bed --keep-allele-order --out ./all_chr.final_rename_${region}.${chr}
        sed -i 's/^/chr/' ./all_chr.final_rename_${region}.${chr}.bim 
    done

    ##### calculate genotype PCs 
    plink --bfile ./all_chr.final_rename_${region} --pca 5 --out ./all_chr.final_rename_${region}.PCs
done 






