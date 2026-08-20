#! /usr/bin/bash
#$ -cwd
#$ -N fdr.v1
#$ -o logv1/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e logv1/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_data=4G,h_rt=2:00:00,highp
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load R

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis/Review_R1/ROSMAP
DIR_CONC=${DIR_WORK}/Concordance
DIR_ASARP=${DIR_WORK}/ASARP/asarp_final
DIR_ASARPOUT=${DIR_WORK}/ASARP/ASARP
events=${DIR_ASARP}/path/to/significant/asAPA_events.txt
DIR_ONE=${DIR_CONC}/STEP1_PREPROC
DIR_TWO=${DIR_CONC}/STEP2_RDD
DIR_FOUR=${DIR_CONC}/STEP4_SISCORES
DIR_FIVE=${DIR_CONC}/STEP5_PEAKS
DIR_SIX=${DIR_CONC}/STEP6_CAUSAL
DIR_SEVEN=${DIR_CONC}/STEP7_NONCAUSAL
DIR_OUT=${DIR_CONC}/STEP8_FDR

###################

echo "--- $(date) --- Correct ---"

method=fdr
minT=minT2

mkdir -p log/fisherP.adjust_r
for indir in ${DIR_SEVEN}/v1
do
	indir1=$indir/Control
	indir2=$indir/AD

	for pthresh in 0.05 0.1
	do
		outdir1=${DIR_OUT}/v1/Control.filtered2.fisherP.$method.$pthresh
		outdir2=${DIR_OUT}/v1/AD.filtered2.fisherP.$method.$pthresh

		outf_Control=$outdir1
		outf_AD=$outdir2

		mkdir -p $outf_Control
                mkdir -p $outf_AD

		R CMD BATCH --no-save '--args indir="'$indir1'" padjmethod="'$method'" pthresh='$pthresh' outdir="'$outf_Control'"' ./scripts/fisherP.adjust.R log/fisherP.adjust_r/.si.detect.non-causal.v1.Control
		R CMD BATCH --no-save '--args indir="'$indir2'" padjmethod="'$method'" pthresh='$pthresh' outdir="'$outf_AD'"' ./scripts/fisherP.adjust.R log/fisherP.adjust_r/.si.detect.non-causal.v1.AD
	done

done

echo "--- $(date) --- Done ---"
