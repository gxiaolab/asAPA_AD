#! /usr/bin/bash
#$ -cwd
#$ -N rosmaprbp1
#$ -o log/job.$JOB_NAME.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_rt=24:00:00,h_data=8G
#$ -pe shared 2
#$ -t 1-1:1

THREAD=2

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python

####################

DIR_BASE=/home/directory/asAPA_analysis
DIR_WORK=${DIR_BASE}/RBP
DIR_OUT=${DIR_WORK}/MAPPED_step1
working=${DIR_WORK}/working2

eCLIPIDs=/path/to/eCLIP/eCLIP_filenames.txt
eclip_dir=/path/to/eCLIP
rbp_script=./rbp_find.py

echo "--- $(date) --- Set Files ---"

out_asapa=${DIR_OUT}/asapa_snps_mapped.txt
out_ctrl=${DIR_OUT}/ctrl_snps_mapped.txt

event_file=${DIR_WORK}/asapa_snps.txt
ctrl_file=${DIR_WORK}/ctrl_snps.txt

echo "--- $(date) --- Parse Control File ---"

python ${rbp_script} \
    ${ctrl_file} \
    ${eCLIPIDs} \
    ${eclip_dir} \
    ${out_ctrl}

echo "--- $(date) --- Parse Event File ---"

python ${rbp_script} \
    ${event_file} \
    ${eCLIPIDs} \
    ${eclip_dir} \
    ${out_asapa}

echo "--- $(date) --- Done ---"
