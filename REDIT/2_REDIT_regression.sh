#! /usr/bin/bash
#$ -cwd
#$ -N reditreg
#$ -o log/job.$JOB_NAME.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_rt=2:00:00,h_data=2G
#$ -pe shared 1
#$ -t 1-1:1 

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load R

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_ASARP=${DIR_BASE}/ASARP
DIR_ASE=${DIR_ASARP}/asarp_final
DIR_PREP=${DIR_OUT}/prep
DIR_STEP1=${DIR_OUT}/STEP1_COUNT
DIR_STEP2=${DIR_OUT}/STEP2_ORGANIZE
DIR_STEP3=${DIR_OUT}/STEP3_REDITTESTABLE

DIR_REVIEW=${DIR_BASE}/asAPA_analysis
DIR_STEP4=${DIR_WORK}/REDIT/STEP4_REDIT
DIR_STEP5=${DIR_WORK}/REDIT/STEP5_PVALCORRECT

events=/path/to/significant_asapa_events.txt
rnaseqids_redit=${DIR_PREP}/rnaseqids.txt
metadata_map=/path/to/metadata.tsv
metadata_cov=/path/to/metadata_covariates.txt
rin_file=/path/to/map/brain_to_RIN.txt
brain_file=/path/to/map/rnaseqid_to_samp.txt
neuron=/path/to/map/all_neuronal_proportions.txt

REDIT=./REDIT_regression.R
fdr_correction=./fdr_correction.py

echo "--- $(date) --- Make Output File ---"

rnaseqid=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${rnaseqids_redit})

echo "--- $(date) --- Step 4: Run REDIT ---"

# Task to 1-# of REDIT testable files

out=${DIR_STEP4}/${rnaseqid}

while read asAPA
do
    rnaid=$(echo -n "${asAPA}" | cut -f 1)
    gene=$(echo -n "${asAPA}" | cut -f 3)
    chrom=$(echo -n "${asAPA}" | cut -f 4)
    pos=$(echo -n "${asAPA}" | cut -f 5)
    tissue=$(echo -n "${asAPA}" | cut -f 15)
    tissue2="${tissue%]}"
    tissue2="${tissue2#[}"
    tissue2=$(echo "${tissue2}" | tr ' ' '_')
    stat=$(echo -n "${asAPA}" | cut -f 17)
    refalt=$(echo -n "${asAPA}" | cut -f 6)
    tmp=${DIR_STEP4}/${rnaid}_${pos}_tmp.txt
    tmp2=${DIR_STEP4}/${rnaid}_${pos}_tmp2.txt

    # Create temporary file
    # Can change count cutoffs! ($7 and $8)
    awk -F'\t' -v gn="${gene}" -v ch="${chrom}" -v ps="${pos}" -v tg="${refalt}" '$3 == gn && $4 == ch && $5 == ps && $6 == tg && $7 >= 0 && $8 >= 0' "${DIR_STEP2}/${tissue2}_AD.txt" >> ${tmp}
    awk -F'\t' -v gn="${gene}" -v ch="${chrom}" -v ps="${pos}" -v tg="${refalt}" '$3 == gn && $4 == ch && $5 == ps && $6 == tg && $7 >= 0 && $8 >= 0' "${DIR_STEP2}/${tissue2}_Control.txt" >> ${tmp}

    # Can change count cutoffs! ($7 and $8)
    ad_count=$(awk -F'\t' -v gn="${gene}" -v ch="${chrom}" -v ps="${pos}" -v tg="${refalt}" '$3 == gn && $4 == ch && $5 == ps && $6 == tg && $7 >= 0 && $8 >= 0' "${DIR_STEP2}/${tissue2}_AD.txt" | wc -l || true)
    ctrl_count=$(awk -F'\t' -v gn="${gene}" -v ch="${chrom}" -v ps="${pos}" -v tg="${refalt}" '$3 == gn && $4 == ch && $5 == ps && $6 == tg && $7 >= 0 && $8 >= 0' "${DIR_STEP2}/${tissue2}_Control.txt" | wc -l || true)
    ad_count_sig=$(awk -F'\t' -v gn="${gene}" -v ch="${chrom}" -v ps="${pos}" -v tg="${refalt}" '$3 == gn && $4 == ch && $5 == ps && $6 == tg && $7 >= 0 && $8 >= 0' "${DIR_STEP2}/${tissue2}_AD.txt" | grep 'Significant' | wc -l || true)
    ctrl_count_sig=$(awk -F'\t' -v gn="${gene}" -v ch="${chrom}" -v ps="${pos}" -v tg="${refalt}" '$3 == gn && $4 == ch && $5 == ps && $6 == tg && $7 >= 0 && $8 >= 0' "${DIR_STEP2}/${tissue2}_Control.txt" | grep 'Significant' | wc -l || true)

    while read entry
    do
        rnaid=$(echo -n "${entry}" | cut -f 1)
        ampad=$(echo -n "${entry}" | cut -f 2)
        gene=$(echo -n "${entry}" | cut -f 3)
        chrom=$(echo -n "${entry}" | cut -f 4)
        pos=$(echo -n "${entry}" | cut -f 5)
        refalt=$(echo -n "${entry}" | cut -f 6)
        ref_count=$(echo -n "${entry}" | cut -f 7)
        alt_count=$(echo -n "${entry}" | cut -f 8)
        tissue=$(echo -n "${entry}" | cut -f 15)
        stat=$(echo -n "${entry}" | cut -f 17)

        age=$(grep "${ampad}" ${metadata_cov} | cut -f 8 | head -n 1)
        if [ "$age" == "90+" ]; then
            age=90
        fi
        sex=$(grep "${ampad}" ${metadata_cov} | cut -f 4 | head -n 1)
        pmi=$(grep "${ampad}" ${metadata_cov} | cut -f 12 | head -n 1)
        race=$(grep "${ampad}" ${metadata_cov} | cut -f 5 | head -n 1)
        yearsEduc=$(grep "${ampad}" ${metadata_cov} | cut -f 7 | head -n 1)
        brainweight=$(grep "${ampad}" ${metadata_cov} | cut -f 14 | head -n 1)

        brain=$(grep "${rnaid}" ${brain_file} | cut -f 2 | head -n 1)
        rin=$(grep -w "${brain}" ${rin_file} | cut -f 2 | head -n 1)
        neuronal_prop=$(grep -w "${brain}" ${neuron} | cut -f 2 | head -n 1)
        echo -e "${rnaid}\t${gene}\t${chrom}\t${pos}\t${refalt}\t${ref_count}\t${alt_count}\t${tissue}\t${stat}\t${age}\t${sex}\t${pmi}\t${rin}\t${neuronal_prop}" >> ${tmp2}
    done < "${tmp}"

    Rscript ${REDIT} "${tmp2}" ${ad_count} ${ctrl_count} ${ad_count_sig} ${ctrl_count_sig} >> ${out}

    rm -f ${tmp}
    rm -f ${tmp2}

done < ${DIR_STEP3}/${rnaseqid}

echo "--- $(date) --- Step 5: Organize ---"

# Task to 1-1

while read id
do
        rnaseqid=$(echo -n "${id}" | cut -d'.' -f1)
        tissue_raw=$(grep "${rnaseqid}" ${events} | cut -f8 | sort -u)
        tissue_1="${tissue_raw#[}"
        tissue_2="${tissue_1%]}"
        tissue_clean="$(echo "${tissue_2}" | tr ' ' '_')"
        if [ "${tissue_clean}" == "frontal_pole" ]; then
                cat ${DIR_STEP4}/${rnaseqid}.txt >> ${DIR_STEP5}/FP_tmp.txt || true
        elif [ "${tissue_clean}" == "inferior_frontal_gyrus" ]; then
                cat ${DIR_STEP4}/${rnaseqid}.txt >> ${DIR_STEP5}/IFG_tmp.txt || true
        elif [ "${tissue_clean}" == "parahippocampal_gyrus" ]; then
                cat ${DIR_STEP4}/${rnaseqid}.txt >> ${DIR_STEP5}/PG_tmp.txt || true
        elif [ "${tissue_clean}" == "superior_temporal_gyrus" ]; then
                cat ${DIR_STEP4}/${rnaseqid}.txt >> ${DIR_STEP5}/STG_tmp.txt || true
        fi
done < ${rnaseqids_redit}

wc -l ${DIR_STEP5}/*_tmp.txt
sort -u ${DIR_STEP5}/FP_tmp.txt -o ${DIR_STEP5}/FP_tmp.txt
sort -u ${DIR_STEP5}/IFG_tmp.txt -o ${DIR_STEP5}/IFG_tmp.txt
sort -u ${DIR_STEP5}/PG_tmp.txt -o ${DIR_STEP5}/PG_tmp.txt
sort -u ${DIR_STEP5}/STG_tmp.txt -o ${DIR_STEP5}/STG_tmp.txt
wc -l ${DIR_STEP5}/*_tmp.txt

echo -e "gene\tchrom\tpos\ttissue\tSNP\ttot_ad\ttot_ctrl\tsig_ad\tsig_ctrl\tconverged_01\tage_p\tpmi_p\trin_p\tneuronal_frac_p\tsex_p\tAD_p\tref_ad\tref_ctrl\talt_ad\talt_ctrl\tfull_pvalues\n$(cat ${DIR_STEP5}/FP_tmp.txt)" > ${DIR_STEP5}/FP.txt
echo -e "gene\tchrom\tpos\ttissue\tSNP\ttot_ad\ttot_ctrl\tsig_ad\tsig_ctrl\tconverged_01\tage_p\tpmi_p\trin_p\tneuronal_frac_p\tsex_p\tAD_p\tref_ad\tref_ctrl\talt_ad\talt_ctrl\tfull_pvalues\n$(cat ${DIR_STEP5}/IFG_tmp.txt)" > ${DIR_STEP5}/IFG.txt
echo -e "gene\tchrom\tpos\ttissue\tSNP\ttot_ad\ttot_ctrl\tsig_ad\tsig_ctrl\tconverged_01\tage_p\tpmi_p\trin_p\tneuronal_frac_p\tsex_p\tAD_p\tref_ad\tref_ctrl\talt_ad\talt_ctrl\tfull_pvalues\n$(cat ${DIR_STEP5}/PG_tmp.txt)" > ${DIR_STEP5}/PG.txt
echo -e "gene\tchrom\tpos\ttissue\tSNP\ttot_ad\ttot_ctrl\tsig_ad\tsig_ctrl\tconverged_01\tage_p\tpmi_p\trin_p\tneuronal_frac_p\tsex_p\tAD_p\tref_ad\tref_ctrl\talt_ad\talt_ctrl\tfull_pvalues\n$(cat ${DIR_STEP5}/STG_tmp.txt)" > ${DIR_STEP5}/STG.txt

rm ${DIR_STEP5}/FP_tmp.txt
rm ${DIR_STEP5}/IFG_tmp.txt
rm ${DIR_STEP5}/PG_tmp.txt
rm ${DIR_STEP5}/STG_tmp.txt

echo "--- $(date) --- Step 6: Pvalue Correction ---"

## change task to 1-1

module load python
names=${DIR_STEP5}/names.txt

while read nm
do
        out=${DIR_STEP5}/${nm}_corrected.txt
        python3 ${fdr_correction} ${DIR_STEP5}/${nm}.txt ${out}
done < ${names}

echo "--- $(date) --- Done ---"
