#! /usr/bin/bash
#$ -cwd
#$ -N r_concrbpmap3
#$ -o log3/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log3/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_rt=1:00:00,highp
#$ -pe shared 1
#$ -t 1-4:1

THREAD=2

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python/3.7.3

####################

DIR_BASE=/home/directory/asAPA_analysis
DIR_WORK=${DIR_BASE}/RBP
DIR_ONE=${DIR_WORK}/MAPPED_step1
DIR_TWO=${DIR_WORK}/ORGANIZE_step2
DIR_OUT=${DIR_WORK}/CORRECT_step3

pval_script=./bootstrap_for_rbp_analysis.py
fdr_correction=./fdr_correction.py

echo "--- $(date) --- Make Content File ---"

tissue="Tissue Type" # Replace with the actual tissue name
out=${DIR_TWO}/${tissue// /_}.txt

echo "--- $(date) --- Get pvalue  ---"

out_2=${DIR_OUT}/${tissue// /_}_withpvalue.txt
plot_path=${DIR_OUT}/example_plots
mkdir -p ${plot_path}
echo -e "RBP\tTest_Count\tTotal_Test\tTest_Prop\tControl_Count\tTotal_Control\tControl_Prop\tProp_Diff\tFold\tPvalue" >> ${out_2}
python3 ${pval_script} ${out} ${plot_path} >> ${out_2}

echo "--- $(date) --- pvalue correction  ---"

out_3=${DIR_OUT}/${tissue// /_}_withpvalue_corrected.txt
python3 ${fdr_correction} ${out_2} ${out_3}
rm ${out_2}

echo "--- $(date) --- Done  ---"
