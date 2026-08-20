#! /usr/bin/bash
#$ -cwd
#$ -N finredit
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=1:00:00,highp
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_OUT=${DIR_WORK}/REDIT
events=/path/to/all_asapa_events.txt
sigevents=/path/to/significant_asapa_events.txt

DIR_PREP=${DIR_OUT}/prep
DIR_STEP1=${DIR_OUT}/STEP1_COUNT
DIR_STEP2=${DIR_OUT}/STEP2_ORGANIZE
DIR_STEP3=${DIR_OUT}/STEP3_REDITTESTABLE

rnaseqids=${DIR_PREP}/rnaseqids.txt

working=${DIR_ALLELE}/working
DIR_RDD=/path/to/rdd/scripts
map=/path/to/map/geneid_to_gene_name.txt

echo "--- $(date) --- Make Output File ---"

rnaseqid=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${rnaseqids})
out1=${DIR_STEP1}/${rnaseqid}.txt

echo "--- $(date) --- Add Counts: Step 1 ---"

# Change task to 1-# of samples with significant asAPA

tmp=${working}/${rnaseqid}_tmp.txt
tmp_unsorted=${working}/${rnaseqid}_tmp_unsorted.txt
grep "${rnaseqid}" ${events} >> ${tmp}

while read event
do
        brain=$(echo -n "${event}" | cut -f 2)
        gene=$(echo -n "${event}" | cut -f 3)
        geneid=$(awk -F'\t' -v gn="${gene}" '$2==gn' ${map} | cut -f 1)
        ase_check=$(grep -w "${geneid}" ${DIR_ASE}/${rnaseqid}.ase.prediction || true)
        if [ ! -z "${ase_check}" ]; then
                echo "ASE GENE"
                continue
        fi

        chrom=$(echo -n "${event}" | cut -f 4)
        pos_targ=$(echo -n "${event}" | cut -f 5)
        nev=$(echo -n "${event}" | cut -f 6)
        sig=$(awk -F'\t' -v rid="${rnaseqid}" -v gn="${gene}" -v ch="${chrom}" -v ps="${pos_targ}" '$1==rid && $3==gn && $4==ch && $5==ps' ${sigevents} || true)
        if [ -z "${sig}" ]; then
                sig="non-significant"
        else
                sig="Significant"
        fi

        tissue=$(echo -n "${event}" | cut -f 8)
        spec_stat=$(echo -n "${event}" | cut -f 9)
        stat=$(echo -n "${event}" | cut -f 10)
        jobid=$(echo -n "${event}" | cut -f 11 | tr -d '\n')
        tmp_targ=${working}/${rnaseqid}_tmp_targsnp.txt
        awk -F' ' -v ch="${chrom}" -v ps="${pos_targ}" '$1==ch && $2==ps {print $3}' ${DIR_RDD}/${rnaseqid}.final.all_snvs.${rnaseqid}.rdd >> ${tmp_targ} || true
        if [ ! -s "${tmp_targ}" ]; then
                echo "NO TAG SNP IN RDD FILE"
        else
                while read tg_snp
                do
                        targ_snp=$(echo -n "${tg_snp}" | cut -f 1)
                        targ_ref=$(awk -F' ' -v ch="${chrom}" -v ps="${pos_targ}" -v snp="${targ_snp}" '$1==ch && $2==ps && $3==snp {print $5}' ${DIR_RDD}/${rnaseqid}.final.all_snvs.${rnaseqid}.rdd | cut -d':' -f1)
                        targ_alt=$(awk -F' ' -v ch="${chrom}" -v ps="${pos_targ}" -v snp="${targ_snp}" '$1==ch && $2==ps && $3==snp {print $5}' ${DIR_RDD}/${rnaseqid}.final.all_snvs.${rnaseqid}.rdd | cut -d':' -f2)
                        echo -e "${rnaseqid}\t${brain}\t${gene}\t${chrom}\t${pos_targ}\t${targ_snp}\t${targ_ref}\t${targ_alt}\tcontrol_pos\tcontrol_refalt\tcontrol_ref_count\tcontrol_alt_count\t${nev}\t${sig}\t${tissue}\t${spec_stat}\t${stat}\t${jobid}" >> ${tmp_unsorted}
                done < ${tmp_targ}
                rm ${tmp_targ}
        fi
done < ${tmp}

sort -u ${tmp_unsorted} >> ${out1}

rm ${tmp}
rm ${tmp_unsorted}

echo "--- $(date) --- Organize into tissue and status: Step 2 ---"

# Change task to 1-1

while read rnaseqid
do
    br=$(grep "${rnaseqid}" ${events} | cut -f 8 | sort -u)
    br="${br%]}"
    br="${br#[}"
    br=$(echo "${br}" | tr ' ' '_')

    st=$(grep "${rnaseqid}" ${events} | cut -f 10 | sort -u)

    out2=${DIR_STEP2}/${br}_${st}.txt

    cat ${DIR_STEP1}/${rnaseqid}.txt >> ${out2}

    sort -u ${out2} -o ${out2}
done < ${rnaseqids}

echo "--- $(date) --- Find Testable Events: Step 3 ---"

## Change task to 1-# of samples with significant asAPA

tmp=${DIR_STEP3}/${rnaseqid}_tmp.txt
grep "${rnaseqid}" ${sigevents} >> ${tmp}

out3=${DIR_STEP3}/${rnaseqid}.txt

while read asAPA
do
    gene=$(echo -n "${asAPA}" | cut -f 3)
    geneid=$(awk -F'\t' -v gn="${gene}" '$2==gn' ${map} | cut -f 1)
    chrom=$(echo -n "${asAPA}" | cut -f 4)
    pos=$(echo -n "${asAPA}" | cut -f 5)
    tissue=$(echo -n "${asAPA}" | cut -f 8)
    tissue2="${tissue%]}"
    tissue2="${tissue2#[}"
    tissue2=$(echo "${tissue2}" | tr ' ' '_')
    stat=$(echo -n "${asAPA}" | cut -f 10)

    refalt=$(echo -n "${asAPA}" | cut -f 13)
    ref_count=$(echo -n "${asAPA}" | cut -f 14 | cut -d':' -f1)
    alt_count=$(echo -n "${asAPA}" | cut -f 14 | cut -d':' -f2)

    # Can change count cutoffs! ($7 and $8)
    ad_count=$(awk -F'\t' -v gn="${gene}" -v ch="${chrom}" -v ps="${pos}" -v tg="${refalt}" '$3 == gn && $4 == ch && $5 == ps && $6 == tg && $7 >= 0 && $8 >= 0' "${DIR_STEP2}/${tissue2}_AD.txt" | wc -l || true)
    ctrl_count=$(awk -F'\t' -v gn="${gene}" -v ch="${chrom}" -v ps="${pos}" -v tg="${refalt}" '$3 == gn && $4 == ch && $5 == ps && $6 == tg && $7 >= 0 && $8 >= 0' "${DIR_STEP2}/${tissue2}_Control.txt" | wc -l || true)

    if [ "${ad_count}" -ge 5 ] && [ "${ctrl_count}" -ge 5 ]; then
        echo -e "${rnaseqid}\tbrainid\t${gene}\t${chrom}\t${pos}\t${refalt}\t${ref_count}\t${alt_count}\tcontrol_pos\tcontrol_refalt\tcontrol_ref_count\tcontrol_alt_count\tNEV\tSignificant\t${tissue}\tspecific_status\t${stat}\tjob" >> ${out3}
    fi

done < ${tmp}

rm ${tmp}

echo "--- $(date) --- Done ---"
