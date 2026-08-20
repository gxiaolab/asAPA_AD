#! /usr/bin/bash
#$ -cwd
#$ -N summarize_candid
#$ -o log/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e log/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_data=4G,h_rt=2:00:00,highp
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
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
DIR_EIGHT=${DIR_CONC}/STEP8_FDR
DIR_OUT=${DIR_CONC}/STEP9_SUMMARIZE

###################

echo "--- $(date) --- Set parameters ---"

n=10
fdr=0.05
out_dir=${DIR_OUT}
mkdir -p $out_dir/AD
mkdir -p $out_dir/Control

v1_indir=${DIR_EIGHT}/v1
v2_indir=${DIR_EIGHT}/v2
v2b_indir=${DIR_EIGHT}/v2b

v1_AD=$v1_indir/AD.filtered2.fisherP.fdr.$fdr
v1_Control=$v1_indir/Control.filtered2.fisherP.fdr.$fdr
v2_AD=$v2_indir/AD.filtered2.fisherP.fdr.$fdr
v2_Control=$v2_indir/Control.filtered2.fisherP.fdr.$fdr
v2b_AD=$v2b_indir/AD.filtered2.fisherP.fdr.$fdr
v2b_Control=$v2b_indir/Control.filtered2.fisherP.fdr.$fdr

tissues=("Tissue Types")

echo "--- $(date) --- Sort unique ---"

# get the uniq candidates and exons
for tissue in "${tissues[@]}"
do
    # AD
    v1_file=$v1_AD/$tissue.min$n.*.txt
    v2_file=$v2_AD/$tissue.min$n.*.txt
    v2b_file=$v2b_AD/$tissue.min$n.*.txt
    outf_AD=$out_dir/AD/$tissue.min$n.txt

    less $v1_file | awk -v OFS='\t' '{print $1,$4,$5}' > $outf_AD
    less $v2_file | awk -v OFS='\t' '{print $1,$4,$5}' >> $outf_AD
    less $v2b_file | awk -v OFS='\t' '{print $1,$4,$5}' >> $outf_AD
    less $outf_AD | awk '{print $1}' | sort -u > $outf_AD.uniq.candid.txt
    less $outf_AD | awk '{print $2}' | sort -u > $outf_AD.uniq.exon.txt

    # Control
    v1_file=$v1_Control/$tissue.min$n.*.txt
    v2_file=$v2_Control/$tissue.min$n.*.txt
    v2b_file=$v2b_Control/$tissue.min$n.*.txt
    outf_Control=$out_dir/Control/$tissue.min$n.txt
    less $v1_file | awk -v OFS='\t' '{print $1,$4,$5}' > $outf_Control
    less $v2_file | awk -v OFS='\t' '{print $1,$4,$5}' >> $outf_Control
    less $v2b_file | awk -v OFS='\t' '{print $1,$4,$5}' >> $outf_Control
    less $outf_Control | awk '{print $1}' | sort -u > $outf_Control.uniq.candid.txt
    less $outf_Control | awk '{print $2}' | sort -u > $outf_Control.uniq.exon.txt
done

echo "--- $(date) --- Done ---"
