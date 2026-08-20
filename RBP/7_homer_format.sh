#! /usr/bin/bash
#$ -cwd
#$ -N homerformat
#$ -o log/job.$JOB_NAME.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_rt=8:00:00,h_data=8G
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_MOT=${DIR_WORK}/HOMER
DIR_OUT=${DIR_MOT}/STEP0_INFO
ctrl_events=/path/to/ctrl_events.txt
asapa_events=/path/to/asapa_events.txt
map=/path/to/geneid/to/gene_name_mapping.txt
gtf=/path/to/gencode.gtf

####################

out_ctrl=${DIR_MOT}/formatting_and_storage/ctrl_events_formatted.txt
out_asapa=${DIR_MOT}/formatting_and_storage/asapa_preferred_events_formatted.txt
out_asapa_alt=${DIR_MOT}/formatting_and_storage/asapa_unpreferred_events_formatted.txt

tmp_ctrl=${DIR_MOT}/formatting_and_storage/tmp_ctrl_events_formatted.txt
tmp_asapa=${DIR_MOT}/formatting_and_storage/tmp_asapa_preferred_events_formatted.txt
tmp_asapa_alt=${DIR_MOT}/formatting_and_storage/tmp_asapa_unpreferred_events_formatted.txt

echo "--- $(date) --- asAPA Ref ---"

while read event
do
        gene=$(echo -n "${event}" | cut -f 3)
        geneid=$(awk -F'\t' -v gn="${gene}" '$2==gn' ${map} | cut -f 1)
        chrom=$(echo -n "${event}" | cut -f 4)
        pos=$(echo -n "${event}" | cut -f 5)
        tissue=$(echo -n "${event}" | cut -f 8)
        statuss=$(echo -n "${event}" | cut -f 10)
        strand=$(awk -F'\t' -v gid="${geneid}" '$11==gid' ${gtf} | cut -f 3 | sort -u)
        ref_alt=$(echo -n "${event}" | cut -f 13)
        ref=$(echo -n "${ref_alt}" | cut -d '>' -f 1)
        alt=$(echo -n "${ref_alt}" | cut -d '>' -f 2)
        counts=$(echo -n "${event}" | cut -f 14)
        ref_count=$(echo -n "${counts}" | cut -d ':' -f 1)
        alt_count=$(echo -n "${counts}" | cut -d ':' -f 2)
        if [[ "${ref_count}" -gt "${alt_count}" ]]; then
                ref_alt="${ref}>${alt}"
                alt_ref="${alt}>${ref}"
        else
                ref_alt="${alt}>${ref}"
                alt_ref="${ref}>${alt}"
        fi
        echo -e "${gene}\t${chrom}:${pos}:${ref_alt}\t${strand}\t${chrom}:${pos}:${pos}\t${tissue}\t${statuss}\t3_prime_UTR_variant" >> ${tmp_asapa}
        echo -e "${gene}\t${chrom}:${pos}:${alt_ref}\t${strand}\t${chrom}:${pos}:${pos}\t${tissue}\t${statuss}\t3_prime_UTR_variant" >> ${tmp_asapa_alt}
done < ${asapa_events}

sort -u ${tmp_asapa} > ${out_asapa}
sort -u ${tmp_asapa_alt} > ${out_asapa_alt}

echo "--- $(date) --- Control ---"

# while read event
# do
#         gene=$(echo -n "${event}" | cut -f 3)
#         chrom=$(echo -n "${event}" | cut -f 1)
#         pos=$(echo -n "${event}" | cut -f 2)
#         ref=$(echo -n "${event}" | cut -f 5)
#         alt=$(echo -n "${event}" | cut -f 6)
#         strand=$(echo -n "${event}" | cut -f 4)
#         tissue=$(echo -n "${event}" | cut -f 8)
#         statuss=$(echo -n "${event}" | cut -f 9)
#         echo -e "${gene}\t${chrom}:${pos}:${ref}>${alt}\t${strand}\t${chrom}:${pos}:${pos}\t${tissue}\t${statuss}\t3_prime_UTR_variant" >> ${tmp_ctrl}
# done < ${ctrl_events}

# sort -u ${tmp_ctrl} > ${out_ctrl}

# #rm -f ${tmp_ctrl} ${tmp_asapa}

echo "--- $(date) --- Done ---"
