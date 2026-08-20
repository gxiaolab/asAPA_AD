#! /usr/bin/bash
#$ -cwd
#$ -N homergetseqs
#$ -o log/job.$JOB_NAME.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_rt=2:00:00,h_data=8G
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
DIR_OUT=${DIR_MOT}/STEP1_GETSEQ
DIR_INFO=${DIR_MOT}/STEP0_INFO

####################

echo "--- $(date) --- Get Sequence Around Variants ---"

window_size=7
mer="8mer_preferred"

for dataset in ${DIR_INFO}/*;
do
	indir=${DIR_INFO}

        if [[ "${dataset}" == *"test"* ]]; then
                out_prefix=${DIR_OUT}/higher_apa_events_${mer}
        elif [[ "${dataset}" == *"ctrl"* ]]; then
                out_prefix=${DIR_OUT}/lower_apa_events_${mer}
        fi

	genome_fasta=/path/to/genome.fa
	python ./get_seqs_around_vars.py -i ${dataset} -g $genome_fasta -w $window_size -o $out_prefix
done

echo "--- $(date) --- Done ---"
