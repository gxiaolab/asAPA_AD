#! /usr/bin/bash
#$ -cwd
#$ -N mot_dens_core
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=12:00:00,h_data=8G
#$ -pe shared 4
#$ -t 1-5:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python
module load bedtools

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_FX=${DIR_WORK}/FXS
DIR_FMRP=${DIR_FX}/FMRP_eCLIP
DIR_MOTIF=${DIR_FX}/FMRP_Motif
DIR_PDB=${DIR_FMRP}/polyA_DB

DIR_FOUR=${DIR_PDB}/4_coverage
DIR_SIX=${DIR_PDB}/6_usage
info=${DIR_FX}/info/sample_to_disease_status.txt
cohorts=${DIR_FOUR}/cohorts.txt

genome_fasta=/path/to/genome/fasta.fa
gene_map=/path/to/map/geneid_to_genename.txt
gene_pred=/path/to/gtf_file
map=${DIR_PDB}/human.PAS.map_withPAS_on_Refseq_NM.txt

echo "--- $(date) --- Helper Function ---"

select_best_region() {
    local -a region_array=("$@")
    declare -A freq_map
    local max_freq=0

    # Count frequencies
    for region in "${region_array[@]}"; do
        if [[ -z "$region" ]] || ! [[ "$region" =~ ^[0-9]+[[:space:]]+[0-9]+$ ]]; then
            continue
        fi
        freq_map["$region"]=$(( ${freq_map["$region"]} + 1 ))
        (( freq_map["$region"] > max_freq )) && max_freq=${freq_map["$region"]}
    done

    # No valid regions
    (( max_freq == 0 )) && { echo ""; return; }

    # Get most common regions
    local -a most_common_regions=()
    for region in "${!freq_map[@]}"; do
        if (( freq_map["$region"] == max_freq )); then
            most_common_regions+=("$region")
        fi
    done

    # If exactly one, return it
    if (( ${#most_common_regions[@]} == 1 )); then
        echo "${most_common_regions[0]}"
        return
    fi

    # Tie-breaker: pick longest
    local longest_region=""
    local max_length=0
    for region in "${most_common_regions[@]}"; do
        local start end length
        start=$(echo "$region" | awk '{print $1}')
        end=$(echo "$region" | awk '{print $2}')
        length=$((end - start))
        if (( length > max_length )); then
            max_length=$length
            longest_region=$region
        fi
    done

    echo "$longest_region"
}

echo "--- $(date) --- Get Core ---"

cohort=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${cohorts})
apa_file=${DIR_SIX}/${cohort}/final_usage_difference_.5.txt

out=${DIR_MOTIF}/density_core/${cohort}_regions_case_utr_txsearch.txt

simple_map=${DIR_MOTIF}/working/${cohort}_simple_map_case_utr_txsearch.txt
less ${map} | tr -d ' ' | cut -f1,2,3,6,9,10 | grep -v "original" >> ${simple_map}

while read apa
do
    chrom=$(echo -n "${apa}" | cut -f1)
    pas_start=$(echo -n "${apa}" | cut -f2)
    if [ "${pas_start}" == "PAS_start" ]; then
        continue
    fi
    pas_end=$(echo -n "${apa}" | cut -f3)
    transcript_start=$(echo -n "${apa}" | cut -f4)
    transcript_end=$(echo -n "${apa}" | cut -f5)
    strand=$(echo -n "${apa}" | cut -f6)
    category=$(echo -n "${apa}" | cut -f14)
    if [ "${strand}" == "+" ]; then
        gene=$(awk -F'\t' -v ch="${chrom}" -v pas="${pas_end}" -v tr="${transcript_end}" -v st="${strand}" '$1==ch && $2==pas && $3==st && $6==tr {print $4}' ${simple_map} | head -n1 || true)
        geneid=$(awk -F'\t' -v gn="${gene}" '$2==gn {print $1}' ${gene_map} | head -n1 || true)
        proximal_pas=${pas_end}
        txend=${transcript_end}
    else
        gene=$(awk -F'\t' -v ch="${chrom}" -v pas="${pas_start}" -v tr="${transcript_start}" -v st="${strand}" '$1==ch && $2==pas && $3==st && $5==tr {print $4}' ${simple_map} | head -n1 || true)
        geneid=$(awk -F'\t' -v gn="${gene}" '$2==gn {print $1}' ${gene_map} | head -n1 || true)
        proximal_pas=${pas_start}
        txend=${transcript_start}
    fi

    tmp=${DIR_MOTIF}/working/${cohort}_tmp_swap2.txt
    if [ "${strand}" == "+" ]; then
        awk -F'\t' -v gid=${geneid} -v txe=${txend} '$5==txe && $11==gid' ${gene_pred} >> ${tmp}
    else
        newtxend=$((txend - 1))
        awk -F'\t' -v gid=${geneid} -v txe=${newtxend} '$4==txe && $11==gid' ${gene_pred} >> ${tmp}
    fi

    # Check if tmp file is empty
    if [ ! -s ${tmp} ]; then
        # If the tmp file is empty, run this backup command
        awk -F'\t' -v gid=${geneid} '$11==gid' ${gene_pred} >> ${tmp}
    fi
    
# ---------------------------------------------- #

    # Initialize
    declare -A pas_exon_region_counts=()
    declare -a pas_exon_region_list=()
    declare region_start=""
    declare region_end=""
    region_found=0

    declare -a backup_caseB1_regions=()
    declare -a backup_caseB2_regions=()

    gencode_strand=""

    while IFS=$'\t' read -r tid tx_chrom tx_strand tx_start tx_end cds_start cds_end exon_count exon_starts_str exon_ends_str gene_id _; do
        IFS=',' read -ra exon_starts <<< "$exon_starts_str"
        IFS=',' read -ra exon_ends <<< "$exon_ends_str"

        gencode_strand="${tx_strand}"

        # For each exon, check if PAS falls within
        for ((i=0; i<${#exon_starts[@]}; i++)); do
            exon_start=${exon_starts[$i]}
            exon_end=${exon_ends[$i]}

            if (( proximal_pas >= exon_start && proximal_pas <= exon_end )); then
                if [ "$strand" == "+" ]; then
                    # Check if CDS ends within this exon
                    if (( cds_end >= exon_start && cds_end <= proximal_pas )); then
                        # 3'UTR starts at CDS end, not exon start
                        utr_start=$cds_end
                    else
                        utr_start=$exon_start
                    fi
                    region="${utr_start} ${proximal_pas}"
                else
                    # Negative strand: CDS start marks 3'UTR boundary
                    if (( cds_start >= proximal_pas && cds_start <= exon_end )); then
                        utr_end=$cds_start
                    else
                        utr_end=$exon_end
                    fi
                    region="${proximal_pas} ${utr_end}"
                fi

                pas_exon_region_counts["$region"]=$(( ${pas_exon_region_counts["$region"]} + 1 ))
                pas_exon_region_list+=("$region")
            fi

        done

        # Case B1: CDS boundary to PAS
        if [ "$strand" == "+" ]; then
            if (( cds_end < proximal_pas )); then
                backup_caseB1_regions+=("${cds_end} ${proximal_pas}")
            fi
        else
            if (( cds_start > proximal_pas )); then
                backup_caseB1_regions+=("${proximal_pas} ${cds_start}")
            fi
        fi

        # Case B2: Closest exon upstream (+) or downstream (-)
        closest_distance=999999999
        closest_region=""

        for ((i=0; i<${#exon_starts[@]}; i++)); do
            exon_start=${exon_starts[$i]}
            exon_end=${exon_ends[$i]}

            if [ "$strand" == "+" ]; then
                if (( exon_end < proximal_pas )); then
                    distance=$((proximal_pas - exon_end))
                    if (( distance < closest_distance )); then
                        closest_distance=$distance
                        closest_region="${exon_start} ${exon_end}"
                    fi
                fi
            else
                if (( exon_start > proximal_pas )); then
                    distance=$((exon_start - proximal_pas))
                    if (( distance < closest_distance )); then
                        closest_distance=$distance
                        closest_region="${exon_start} ${exon_end}"
                    fi
                fi
            fi
        done

        if [[ -n "$closest_region" ]]; then
            backup_caseB2_regions+=("$closest_region")
        fi

    done < ${tmp}

# --- REGION SELECTION LOGIC --- #

    # Case A: PAS falls within an exon in at least one transcript
        # Case A: PAS falls within an exon in at least one transcript
    if (( region_found == 0 )); then
    if (( ${#pas_exon_region_list[@]} > 0 )); then
        # Step 1: Find the max frequency
        max_count=0
        for region in "${!pas_exon_region_counts[@]}"; do
            count=${pas_exon_region_counts[$region]}
            (( count > max_count )) && max_count=$count
        done

        # Step 2: Collect all tied regions
        declare -a tied_regions=()
        for region in "${!pas_exon_region_counts[@]}"; do
            if (( pas_exon_region_counts[$region] == max_count )); then
                tied_regions+=("$region")
            fi
        done

        # Step 3: Resolve ties
        if (( ${#tied_regions[@]} == 1 )); then
            selected_region=${tied_regions[0]}
        else
            # Try to compute intersection across tied regions
            region_start=$(printf "%s\n" "${tied_regions[@]}" | awk '{print $1}' | sort -n | tail -1)
            region_end=$(printf "%s\n" "${tied_regions[@]}" | awk '{print $2}' | sort -n | head -1)

            if (( region_start < region_end )); then
                # Non-empty intersection → take it
                selected_region="$region_start $region_end"
            else
                # Empty intersection → fall back to shortest region
                shortest_len=999999999
                selected_region=""
                for region in "${tied_regions[@]}"; do
                    start=$(echo $region | awk '{print $1}')
                    end=$(echo $region | awk '{print $2}')
                    len=$((end - start))
                    if (( len < shortest_len )); then
                        shortest_len=$len
                        selected_region=$region
                    fi
                done
            fi
        fi

        # Step 4: Assign final coordinates
        region_found=1
        region_start=$(echo $selected_region | awk '{print $1}')
        region_end=$(echo $selected_region | awk '{print $2}')
        case="Case_1"
    fi
    fi


    # Case B1: coding region to PAS
        # Case B1: coding region to PAS
    if (( region_found == 0 && ${#backup_caseB1_regions[@]} > 0 )); then
        # Count region frequencies
        declare -A region_counts=()
        max_count=0
        for region in "${backup_caseB1_regions[@]}"; do
            region_counts["$region"]=$(( ${region_counts["$region"]} + 1 ))
            (( region_counts["$region"] > max_count )) && max_count=${region_counts["$region"]}
        done

        # Collect tied regions
        declare -a tied_regions=()
        for region in "${!region_counts[@]}"; do
            if (( region_counts[$region] == max_count )); then
                tied_regions+=("$region")
            fi
        done

        # Resolve ties
        if (( ${#tied_regions[@]} == 1 )); then
            selected_region=${tied_regions[0]}
        else
            region_start=$(printf "%s\n" "${tied_regions[@]}" | awk '{print $1}' | sort -n | tail -1)
            region_end=$(printf "%s\n" "${tied_regions[@]}" | awk '{print $2}' | sort -n | head -1)

            if (( region_start < region_end )); then
                selected_region="$region_start $region_end"
            else
                shortest_len=999999999
                selected_region=""
                for region in "${tied_regions[@]}"; do
                    start=$(echo $region | awk '{print $1}')
                    end=$(echo $region | awk '{print $2}')
                    len=$((end - start))
                    if (( len < shortest_len )); then
                        shortest_len=$len
                        selected_region=$region
                    fi
                done
            fi
        fi

        region_start=$(echo $selected_region | awk '{print $1}')
        region_end=$(echo $selected_region | awk '{print $2}')
        region_found=1
        case="Case_2"
    fi

    # Case B2: closest exon upstream/downstream
        # Case B2: closest exon upstream (+) or downstream (-)
    if (( region_found == 0 && ${#backup_caseB2_regions[@]} > 0 )); then
        # Count region frequencies
        declare -A region_counts=()
        max_count=0
        for region in "${backup_caseB2_regions[@]}"; do
            region_counts["$region"]=$(( ${region_counts["$region"]} + 1 ))
            (( region_counts["$region"] > max_count )) && max_count=${region_counts["$region"]}
        done

        # Collect tied regions
        declare -a tied_regions=()
        for region in "${!region_counts[@]}"; do
            if (( region_counts[$region] == max_count )); then
                tied_regions+=("$region")
            fi
        done

        # Resolve ties
        if (( ${#tied_regions[@]} == 1 )); then
            selected_region=${tied_regions[0]}
        else
            region_start=$(printf "%s\n" "${tied_regions[@]}" | awk '{print $1}' | sort -n | tail -1)
            region_end=$(printf "%s\n" "${tied_regions[@]}" | awk '{print $2}' | sort -n | head -1)

            if (( region_start < region_end )); then
                selected_region="$region_start $region_end"
            else
                shortest_len=999999999
                selected_region=""
                for region in "${tied_regions[@]}"; do
                    start=$(echo $region | awk '{print $1}')
                    end=$(echo $region | awk '{print $2}')
                    len=$((end - start))
                    if (( len < shortest_len )); then
                        shortest_len=$len
                        selected_region=$region
                    fi
                done
            fi
        fi

        region_start=$(echo $selected_region | awk '{print $1}')
        region_end=$(echo $selected_region | awk '{print $2}')
        region_found=1
        case="Case_3"
    fi


    # Output if valid region found
    if (( region_found == 1 )); then
        echo -e "${chrom}\t${gene}\t${proximal_pas}\t${region_start}\t${region_end}\t${txend}\t${category}\t${strand}\t${gencode_strand}\t${case}" >> ${out}
    fi

# ------------------------------------------- #

    rm ${tmp}

done < ${apa_file}

rm ${simple_map}

echo "--- $(date) --- Done ---"
