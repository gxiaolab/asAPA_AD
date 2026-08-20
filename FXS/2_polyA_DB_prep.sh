#! /usr/bin/bash
#$ -cwd
#$ -N polyadb_prep
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=4:00:00,h_data=4G,highp
#$ -pe shared 1
#$ -t 1-22:1

THREAD=2

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_FX=${DIR_WORK}/FXS
DIR_FMRP=${DIR_FX}/FMRP_eCLIP
DIR_PDB=${DIR_FMRP}/polyA_DB

prep=./polyA_DB_prep.py

echo "--- $(date) --- Expand Positions ---"

in=${DIR_PDB}/2_3utr_by_chrom/chr${idx}_3utr_annotated_human.PAS_on_Refseq.txt
pas=${DIR_PDB}/3_region_pas/chr${idx}_region_annotated_human.PAS_on_Refseq.txt
endt=${DIR_PDB}/3_region_endt/chr${idx}_region_annotated_human.PAS_on_Refseq.txt

python ${prep} ${in} ${pas} ${endt}

echo "--- $(date) --- Done ---"
