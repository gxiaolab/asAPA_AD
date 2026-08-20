#! /usr/bin/bash
#$ -cwd
#$ -N deepripeF
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=4:00:00,highp
#$ -pe shared 1
#$ -t 1-1:1

THREAD=2

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load htslib

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
ctrl_events=${DIR_WORK}/unannotated_ctrl_events.txt
vcf_files=${DIR_WORK}/VCF
vcf_names=${vcf_files}/vcf_names.txt
DIR_LOG=${DIR_WORK}/ASARP_FINAL/log
working=${DIR_WORK}/DeepRiPe/working

echo "--- $(date) --- Make Content File ---"

out=${DIR_WORK}/DeepRiPe/variants/ctrl_events.txt

echo "--- $(date) --- Label Info ---"

while read event
do
        rnaseqid=$(echo -n "${event}" | cut -f 1)
        brain=$(echo -n "${event}" | cut -f 2)
        gene=$(echo -n "${event}" | cut -f 3)
        chrom=$(echo -n "${event}" | cut -f 4)
        pos=$(echo -n "${event}" | cut -f 5)
        nev=$(echo -n "${event}" | cut -f 6)
        tissue=$(echo -n "${event}" | cut -f 8)
        statuss=$(echo -n "${event}" | cut -f 10)
        job=$(echo -n "${event}" | cut -f 11)
        chrom_number=$(echo ${chrom:3})
        if [ "${chrom_number}" == "23" ]
        then
                continue
        fi
        #####################################
        tmp=${working}/${chrom}.${pos}.working.txt
        tmp2=${working}/${chrom}.${pos}.working2.txt
        grep "TESTABLE: " ${DIR_LOG}/${job} | sort -u >> ${tmp}
        awk '!seen[$8, $9]++' ${tmp} >> ${tmp2}
	jobline=$(awk -v chr=${chrom} -v ps=${pos} '$9==chr && $8==ps {print $0}' ${tmp2})
        strand=$(echo "$jobline" | grep -o 'AT: 3[+-][(][0-9.]*[)]' | awk -F'[0-9.]*' '{print $2}' | awk -F'[(]' '{print $1}' || true)
        if [ -z "${strand}" ]
        then
                strand=$(echo "$jobline" | grep -o '[$][+-]>[)][)]' | awk -F'[$]' '{print $2}' | awk -F'>' '{print $1}' || true)
        fi
        rm ${tmp}
        rm ${tmp2}
        #####################################
        corresponding_vcf=${vcf_files}/${chrom_number}.vcf.gz
        search=$(echo "${chrom}:${pos}-${pos}")
        info=$(tabix ${corresponding_vcf} ${search})
        if [ -z "${info}" ]
        then
                echo -e "${chrom}\t${pos}\t${gene}\t${nev}\t${strand}\t${rnaseqid}\t${brain}\t${tissue}\t${statuss}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> ${out}
        else
                reference=$(echo ${info} | cut -d' ' -f 4)
                alternate=$(echo ${info} | cut -d' ' -f 5)
                annot=$(echo ${info} | cut -d' ' -f 8)
                list_kinds=$(echo ${annot} | cut -d'|' -f 2)
                fields=$(echo "${annot}" | grep -oP "(?<=,).*?(?=(,|$))" || true)
                list_more=$(echo "${fields}" | cut -d'|' -f 2)
                list_kinds+=" ${list_more}"
                gene=$(echo ${annot} | cut -d'|' -f 4)
                nmd="NA"
                ## Find kind of annotation first ## 
                count=$(echo "${list_kinds}" | wc -w)
                if [[ $count -eq 1 ]]; then
                        kind=$(echo "$list_kinds" | xargs echo -n)
                elif [[ $count -eq 0 ]]; then
                        kind="NA"
                else
                        if [[ ${list_kinds} == *"3_prime_UTR_variant"* ]]; then
                                kind="3_prime_UTR_variant"
                        elif [[ ${list_kinds} == *"missense_variant&splice_region_variant"* ]]; then
                                kind="missense_variant&splice_region_variant"
                        elif [[ ${list_kinds} == *"missense_variant"* ]]; then
                                kind="missense_variant"
                        elif [[ ${list_kinds} == *"splice_region_variant&synonymous_variant"* ]]; then
                                kind="splice_region_variant&synonymous_variant"
                        elif [[ ${list_kinds} == *"stop_gained"* ]]; then
                                kind="stop_gained"
                        elif [[ ${list_kinds} == *"stop_retained_variant"* ]]; then
                                kind="stop_retained_variant"
                        elif [[ ${list_kinds} == *"synonymous_variant"* ]]; then
                                kind="synonymous_variant"
                        elif [[ ${list_kinds} == *"splice_region_variant&intron_variant"* ]]; then
                                kind="splice_region_variant&intron_variant"
                        elif [[ ${list_kinds} == *"splice_region_variant&non_coding_exon_variant"* ]]; then
                                kind="splice_region_variant&non_coding_exon_variant"
                        elif [[ ${list_kinds} == *"splice_region_variant"* ]]; then
                                kind="splice_region_variant"
                        elif [[ ${list_kinds} == *"splice_donor_variant&intron_variant"* ]]; then
                                kind="splice_donor_variant&intron_variant"
                        elif [[ ${list_kinds} == *"intron_variant"* ]]; then
                                kind="intron_variant"
                        elif [[ ${list_kinds} == *"5_prime_UTR_premature_start_codon_gain_variant"* ]]; then
                                kind="5_prime_UTR_premature_start_codon_gain_variant"
                        elif [[ ${list_kinds} == *"5_prime_UTR_variant"* ]]; then
                                kind="5_prime_UTR_variant"
                        elif [[ ${list_kinds} == *"downstream_gene_variant"* ]]; then
                                kind="downstream_gene_variant"
                        elif [[ ${list_kinds} == *"upstream_gene_variant"* ]]; then
                                kind="upstream_gene_variant"
                        elif [[ ${list_kinds} == *"intergenic_region"* ]]; then
                                kind="intergenic_region"
                        elif [[ ${list_kinds} == *"intragenic_variant"* ]]; then
                                kind="intragenic_variant"
                        elif [[ ${list_kinds} == *"protein_protein_contact"* ]]; then
                                kind="protein_protein_contact"
                        elif [[ ${list_kinds} == *"non_coding_exon_variant"* ]]; then
                                kind="non_coding_exon_variant"
                        else
                                kind="NA"
                        fi
                fi
                ## CADD next ##
                if [[ "$annot" =~ CADD_phred=([0-9]+(\.[0-9]+)?) ]]; then
                        cadd="${BASH_REMATCH[1]}"
                else
                        cadd="NA"
                fi
                ## CLNSIG next ##
                if [[ "$annot" =~ CLNSIG=([0-9]+) ]]; then
                        clnsig="${BASH_REMATCH[1]}"
                else
                        clnsig="NA"
                fi
                ## overlap next ##
                if [[ "$annot" =~ SINE ]]; then
                        overlap="SINE"
                elif [[ "$annot" =~ LINE ]]; then
                        overlap="LINE"
                elif [[ "$annot" =~ LTR ]]; then
                        overlap="LTR"
                elif [[ "$annot" =~ Simple_repeat ]]; then
                        overlap="Simple_repeat"
                else
                        overlap="NA"
                fi
                echo -e "${chrom}\t${pos}\t${gene}\t${nev}\t${strand}\t${rnaseqid}\t${brain}\t${tissue}\t${statuss}\t${reference}\t${alternate}\t${kind}\t${gene}\t${cadd}\t${clnsig}\t${overlap}\t${nmd}" >> ${out}

        fi
done < ${ctrl_events}

echo "--- $(date) --- Done ---"
