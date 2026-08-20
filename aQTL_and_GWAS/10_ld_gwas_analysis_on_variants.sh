#! /usr/bin/bash
#$ -cwd
#$ -N ldgwas
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=12:00:00,h_data=8G
#$ -pe shared 1
#$ -t 1-11:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_GWAS=${DIR_WORK}/GWAS
DIR_STEP0=${DIR_GWAS}/STEP0_CLEANINFO
DIR_STEP1=${DIR_GWAS}/STEP1_LD
DIR_STEP2=${DIR_GWAS}/STEP2_LD200kb
DIR_STEP3=${DIR_GWAS}/STEP3_LDexpand200kb
DIR_STEP4=${DIR_GWAS}/STEP4_WHOLECATALOG
DIR_STEP5=${DIR_GWAS}/STEP5_WHOLECATALOG_ALLCTRL
DIR_STEP6=${DIR_GWAS}/STEP6_PVALUE
DIR_INFO=${DIR_GWAS}/My_SNPs
DIR_OG=${DIR_GWAS}/GWAS_SNPs

# Change
events=${DIR_INFO}/asapa_snps_corrected.txt
event_prefix="asapa"
gwas_filenames=${DIR_STEP0}/gwas_filenames.txt  # List of download GWAS SNP information files
ld_script=./ld_get_variants_in.py

asapa=${DIR_BASE}/path/to/all_significant/asapa_events.txt

echo "--- $(date) --- Get LD info ---"

gwas_file=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${gwas_filenames})

disease="${gwas_file%.txt}"
gwas_events=${DIR_STEP0}/${gwas_file}

# Change
out=${DIR_STEP6}/${disease}_and_${event_prefix}_overlap.txt
working=${DIR_STEP6}/working_${disease}_${event_prefix}
mkdir -p ${working}

while read event; do
    chrom=$(echo -n "${event}" | cut -f 1)
    position=$(echo -n "${event}" | cut -f 2 | tr -d $'\n')

    gene=$(awk -F'\t' -v ch="${chrom}" -v ps="${position}" '$4==ch && $5==ps {print $3}' ${asapa} | sort -u)

    tmp=${working}/${chrom}_${position}_tmp.txt
    python3 ${ld_script} "${event}" >> ${tmp}
    sort -o ${tmp} ${tmp}

    # Find overlapping GWAS SNPs
    matched=$(comm -12 ${tmp} ${gwas_events})

    # For pval
    exact=$(awk -F'\t' -v ch="${chrom}" -v ps="${position}" '$1==ch && $2==ps {print "Exact"}' ${gwas_events})
    if [ -z "${exact}" ]; then
        exact="NA"
    fi

    match_found=0  # Flag to track success

    # Loop through each overlapping GWAS SNP
    while read gwas_snp; do
        gwas_chr=$(echo "$gwas_snp" | cut -f1)
        gwas_pos=$(echo "$gwas_snp" | cut -f2)

        # Check if same chromosome and within 200kb of original event SNP
        if [[ "$gwas_chr" == "$chrom" ]]; then
            diff=$(( position > gwas_pos ? position - gwas_pos : gwas_pos - position ))
            if [ "$diff" -le 200000 ]; then
                ch_num="${gwas_chr#chr}"
                search="${ch_num}:${gwas_pos}"
                pval=$(grep "${search}" ${DIR_OG}/${disease}.tsv | cut -f2)
                echo -e "${chrom}\t${position}\t${gwas_chr}\t${gwas_pos}\t${gene}\t${pval}\t${exact}" >> ${out}
            fi
        fi
    done <<< "$matched"
    rm ${tmp}
done < ${events}

rm -r ${working}

echo "--- $(date) --- Done ---"
