#! /usr/bin/bash
#$ -cwd
#$ -N norm_GO
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=1:00:00,h_data=4G
#$ -pe shared 2
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load R

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_GO=${DIR_WORK}/GO

go_file=${DIR_GO}/genes_of_interest.txt
go=./perform_GO_enrichment.gene_length_RPKM_controlled.normal_distribution.R

echo "--- $(date) --- Normalize ---"

Rscript ${go} ${go_file} /path/to/TPM/average_tpm.txt

echo "--- $(date) --- Done ---"
