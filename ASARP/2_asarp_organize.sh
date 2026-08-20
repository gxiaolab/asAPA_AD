#! /usr/bin/bash
#$ -cwd
#$ -N asarp_organize
#$ -o log/job.JOBNAME.TASK_ID.out
#$ -e log/job.JOBNAME.TASK_ID.err
#$ -V
#$ -l h_rt=2:00:00,h_data=2G,highp
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -e -x -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python
module load R

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/ASARP
DIR_FQ=${DIR_WORK}/FQ
DIR_BAM=${DIR_WORK}/BAM
DIR_PSI=${DIR_WORK}/PSI
DIR_RDD=${DIR_WORK}/RDD
DIR_VCF=${DIR_RDD}/VCF
DIR_CONFIG=${DIR_WORK}/CONFIG
DIR_ASARP=${DIR_WORK}/ASARP
DIR_HIST=${DIR_WORK}/HIST

ALLFQS=${DIR_FQ}/samples_names.txt

map=/map/individual/to/vcf
genemap=/map/genename/to/geneid
metadata=/path/to/metadata/file
metadata_clinical=/path/to/metadata/clinical/file
fqfilename=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${ALLFQS})
rnaseqid=$(awk -v lin=${idx} 'FNR == (lin) {print $2;}' ${ALLFQS})
indiv=$(awk -v lin=${idx} 'FNR == (lin) {print $3;}' ${ALLFQS})
vcf_sample_id=$(grep -w "${indiv}" ${map} | cut -f2 | head -n1)

DIR_OUT=${DIR_WORK}/asarp_final
logfile=${DIR_WORK}/logasarp/job.asarp.${idx}.out
asarp_file=${DIR_ASARP}/${rnaseqid}.out_asarp
rdd_file=${DIR_RDD}/${rnaseqid}.final.all_snvs.${rnaseqid}.rdd

out_all=${DIR_OUT}/${rnaseqid}_all.txt
out_sig=${DIR_OUT}/${rnaseqid}_sig_labeled.txt
out_sig_all=${DIR_OUT}/${rnaseqid}_sig_labeled_all.txt

echo "--- $(date) --- Organize all output  ---"

tmp1=${DIR_OUT}/${rnaseqid}_tmp1.txt
tmp2=${DIR_OUT}/${rnaseqid}_tmp2.txt
tmp3=${DIR_OUT}/${rnaseqid}_tmp3.txt

rm -f ${tmp1} ${tmp2} ${tmp3} ${out_all} ${out_sig} ${out_sig_all}

awk '
    /^Alternative Termination$/ {in_term=1; next}
    in_term && ($0 ~ /^Alternative / || 0~/GENES/ || 0~/SNVs/ || $0 ~ /^There are / || $0 ~ /^#/) {in_term=0}
    in_term && $0 ~ /^ENSG/ {print}
' "${asarp_file}" > "${tmp1}"

tissue="[Tissue_Type(s)]"
log="job.asarp.${idx}.out"

msex=$(awk -F'\t' -v iid="${indiv}" -v syn="${rnaseqid}" '$3 == iid && $2 == syn {print $6}' ${metadata})
    if [ "${msex}" == "1" ]; then
        sex="male"
    elif [ "${msex}" == "0" ]; then
        sex="female"
    else
        sex="NA"
    fi

AoD=$(awk -F',' -v iid="${indiv}" '$18 == iid {print $10}' "${metadata_clinical}")
pmi=$(awk -F',' -v iid="${indiv}" '$18 == iid {print $13}' "${metadata_clinical}")
Braak=$(awk -F',' -v iid="${indiv}" '$18 == iid {print $14}' "${metadata_clinical}")
CERAD=$(awk -F',' -v iid="${indiv}" '$18 == iid {print $15}' "${metadata_clinical}")
Cogdx=$(awk -F',' -v iid="${indiv}" '$18 == iid {print $16}' "${metadata_clinical}")

status=$(Rscript ./organize2.R "${Braak}" "${CERAD}" "${Cogdx}" 2>/dev/null | tail -n1 | tr -d '\r' | xargs)
case "${status}" in
    Not|Low)
        diagnosis="Control"
        ;;
    Intermediate|High)
        diagnosis="AD"
        ;;
    *)
        diagnosis="undetermined"
        ;;
esac

grep "TESTABLE: " ${logfile} | grep "AT: 3" >> ${tmp2}

python ./organize1.py ${tmp2} ${tmp3}

while read logevent
do
    current_signif=$(echo -n "${logevent}" | cut -f1 | cut -d' ' -f1)
    if [ "${current_signif}" == "signif" ]; then
        continue
    fi
    gene_id=$(echo -n "${logevent}" | cut -f2)
    gene=$(awk -F'\t' -v gid="${gene_id}" '$1 == gid {print $2}' ${genemap})
    if [ -z "${gene}" ]; then
        gene="NA"
    fi
    pos=$(echo -n "${logevent}" | cut -f3)
    chrom=$(echo -n "${logevent}" | cut -f4)
    at_sign=$(echo -n "${logevent}" | cut -f5)
    NEV=$(echo -n "${logevent}" | cut -f6)
    ctrl_pos=$(echo -n "${logevent}" | cut -f7)
    SNP=$(echo -n "${logevent}" | cut -f8)
    ref_ctrl=$(echo -n "${logevent}" | cut -f9)
    alt_ctrl=$(echo -n "${logevent}" | cut -f10)
    ref=$(awk -v chr="${chrom}" -v pos="${pos}" '$1 == chr && $2 == pos {print $5; exit}' ${rdd_file} | cut -d':' -f1)
    alt=$(awk -v chr="${chrom}" -v pos="${pos}" '$1 == chr && $2 == pos {print $5; exit}' ${rdd_file} | cut -d':' -f2)
    echo -e "${rnaseqid}\t${indiv}\t${gene}\t${chrom}\t${pos}\t${NEV}\t${current_signif}\t${tissue}\t${status}\t${diagnosis}\t${sex}\t${log}\t${SNP}\t${ref}:${alt}\t${gene_id}\t${ctrl_pos}\t${ref_ctrl}\t${alt_ctrl}\t${AoD}\t${pmi}\t${Braak}\t${CERAD}\t${Cogdx}\t${at_sign}" >> ${out_all}
done < ${tmp3}

echo "--- $(date) --- Organize Sig output  ---"

while read event
do
    rnaseqid=$(echo -n "${event}" | cut -f1)
    indiv=$(echo -n "${event}" | cut -f2)
    gene=$(echo -n "${event}" | cut -f3)
    chrom=$(echo -n "${event}" | cut -f4)
    pos=$(echo -n "${event}" | cut -f5)
    NEV=$(echo -n "${event}" | cut -f6)
    signif=$(echo -n "${event}" | cut -f7)
    tissue=$(echo -n "${event}" | cut -f8)
    status=$(echo -n "${event}" | cut -f9)
    diagnosis=$(echo -n "${event}" | cut -f10)
    sex=$(echo -n "${event}" | cut -f11)
    log=$(echo -n "${event}" | cut -f12)
    SNP=$(echo -n "${event}" | cut -f13)
    ref_alt=$(echo -n "${event}" | cut -f14)
    geneid=$(echo -n "${event}" | cut -f15)
    ctrl_pos=$(echo -n "${event}" | cut -f16)
    ref_ctrl=$(echo -n "${event}" | cut -f17)
    alt_ctrl=$(echo -n "${event}" | cut -f18)
    AoD=$(echo -n "${event}" | cut -f19)
    pmi=$(echo -n "${event}" | cut -f20)
    Braak=$(echo -n "${event}" | cut -f21)
    CERAD=$(echo -n "${event}" | cut -f22)
    Cogdx=$(echo -n "${event}" | cut -f23)
    strand=$(echo -n "${event}" | cut -f24)
   
    new_signif="non-significant"
    if [[ "${ref_ctrl}" =~ ^[0-9]+$ ]] && [[ "${alt_ctrl}" =~ ^[0-9]+$ ]]; then
        total_ctrl=$((ref_ctrl + alt_ctrl))


        if [ "${total_ctrl}" -gt 0 ]; then
            ctrl_ratio=$(awk -v r="${ref_ctrl}" -v t="${total_ctrl}" 'BEGIN{print r/t}')


            ratio_ok=$(awk -v x="${ctrl_ratio}" 'BEGIN{if (x >= 0.45 && x <= 0.55) print "yes"; else print "no"}')


            nev_ok=$(awk -v x="${NEV}" 'BEGIN{if (x <= 0.8) print "yes"; else print "no"}')


            if [ "${ratio_ok}" == "yes" ] && [ "${nev_ok}" == "yes" ]; then
                match_found=$(awk -v gid="${geneid}" -v chr="${chrom}" -v p="${pos}" -v snp="${SNP}" '
                    $1 == gid && $2 == chr && $3 == p && $5 == snp {print "yes"; exit}
                ' "${tmp1}")


                if [ "${match_found}" == "yes" ]; then
                    new_signif="Significant"


                    ra=$(awk -v gid="${geneid}" -v chr="${chrom}" -v p="${pos}" -v snp="${SNP}" '
                        $1 == gid && $2 == chr && $3 == p && $5 == snp {print $6; exit}
                    ' "${tmp1}")


                    ref=$(echo -n "${ra}" | cut -d':' -f1)
                    alt=$(echo -n "${ra}" | cut -d':' -f2)
                    ref_alt="${ref}:${alt}"
                fi
            fi
        fi
    fi
    echo -e "${rnaseqid}\t${indiv}\t${gene}\t${chrom}\t${pos}\t${NEV}\t${new_signif}\t${tissue}\t${status}\t${diagnosis}\t${sex}\t${log}\t${SNP}\t${ref_alt}\t${gene_id}\t${ctrl_pos}\t${ref_ctrl}\t${alt_ctrl}\t${AoD}\t${pmi}\t${Braak}\t${CERAD}\t${Cogdx}\t${strand}" >> ${out_sig_all}
    echo -e "${rnaseqid}\t${indiv}\t${gene}\t${chrom}\t${pos}\t${NEV}\t${new_signif}\t${tissue}\t${status}\t${diagnosis}\t${sex}\t${log}\t${SNP}\t${ref_alt}" >> ${out_sig}
done < ${out_all}

rm -f ${tmp1} ${tmp2} ${tmp3}

echo "--- $(date) --- Done ---"