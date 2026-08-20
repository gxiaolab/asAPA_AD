#! /usr/bin/bash
#$ -cwd
#$ -N deepripeS
#$ -o log3/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log3/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=8:00:00,h_data=16G
#$ -pe shared 4
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load bedtools

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_OUT=${DIR_WORK}/RBP/DeepRiPe

DIR_END=${DIR_OUT}/score_distributions
mkdir -p ${DIR_END}

DIR_TOOLS=$path/to/DeepRiPe/scripts
genome=/path/to/genome/hg38.fa
plot_dir=${DIR_OUT}/plots
score=./deepripe_score_variants.py

inf=${DIR_OUT}/variants/snps.bed
models=${DIR_OUT}/models/deepripe_models.txt

####################

## For CTRL files
names=${DIR_OUT}/variants/ctrl_samples.txt
bedname=$(awk -v lin=${idx} 'FNR == (lin) {print $1;}' ${names})

echo "--- $(date) --- Run ---"

for tr in 0.0;
do
        outf=${DIR_END}/withpref_deepripe_scored_variants.txt	## FOR TEST
        outf=${DIR_END}/ctrl${idx}_deepripe_scored_variants.txt	## FOR CONTROLS
        python ${score} -v $inf -m $models -g $genome -p $plot_dir -t $tr -o $outf
done

echo "--- $(date) --- Done ---"
