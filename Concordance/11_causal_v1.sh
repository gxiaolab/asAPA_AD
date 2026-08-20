#! /usr/bin/bash
#$ -cwd
#$ -N causal_v1
#$ -o logv1/job.$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e logv1/job.$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -V
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load python/2.7.18

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
DIR_OUT=${DIR_CONC}/STEP6_CAUSAL/v1

###################

echo "--- $(date) --- Set variables ---"

mkdir -p ${DIR_OUT}/AD
mkdir -p ${DIR_OUT}/Control

tissue="Tissue Type"

tissue_Control=${tissue}/Control
tissue_AD=${tissue}/AD
annoI=/path/to/gencode.v38.GRCh38.intron.bed
annoE=/path/to/gencode.v38.GRCh38.exon.bed
indir=${DIR_FIVE}/peaks.gmm
ref=/path/to/analysis/Concordance/STEP5_PEAKS/files/apa_r.other.gt.txt
outf_Control=${DIR_OUT}/Control/${tissue}
outf_AD=${DIR_OUT}/AD/${tissue}
si=0.8
pval=0.1
major=0.9

si_indir=${DIR_FOUR}/si.results

echo "--- $(date) --- Get Causal ---"

        si_Controln=$(ls $si_indir/${tissue}/Control/chr3* | wc -l)
        si_ADn=$(ls $si_indir/${tissue}/AD/chr3* | wc -l)
        
        minPt5Control=$(bc <<< 0.05*$si_Controln/1)
        minPt5AD=$(bc <<< 0.05*$si_ADn/1)
        
        minPt10Control=$(bc <<< 0.1*$si_Controln/1)
        minPt10AD=$(bc <<< 0.1*$si_ADn/1)
        
        minPt20Control=$(bc <<< 0.2*$si_Controln/1)
        minPt20AD=$(bc <<< 0.2*$si_ADn/1)

        python2 ./get.causal.v1.py -i $annoI -e $annoE -r $indir -t $tissue_Control -p $pval -n $minPt5Control -m $major -s $si -o $outf_Control.min5.$minPt5Control.txt  
        python2 ./get.causal.v1.py -i $annoI -e $annoE -r $indir -t $tissue_AD -p $pval -n $minPt5AD -m $major -s $si -o $outf_AD.min5.$minPt5AD.txt

        python2 ./get.causal.v1.py -i $annoI -e $annoE -r $indir -t $tissue_Control -p $pval -n $minPt10Control -m $major -s $si -o $outf_Control.min10.$minPt10Control.txt  
        python2 ./get.causal.v1.py -i $annoI -e $annoE -r $indir -t $tissue_AD -p $pval -n $minPt10AD -m $major -s $si -o $outf_AD.min10.$minPt10AD.txt

        python2 ./get.causal.v1.py -i $annoI -e $annoE -r $indir -t $tissue_Control -p $pval -n $minPt20Control -m $major -s $si -o $outf_Control.min20.$minPt20Control.txt  
        python2 ./get.causal.v1.py -i $annoI -e $annoE -r $indir -t $tissue_AD -p $pval -n $minPt20AD -m $major -s $si -o $outf_AD.min20.$minPt20AD.txt

        python2 ./get.causal.v1.py -i $annoI -e $annoE -r $indir -t $tissue -p $pval -n 10 -m $major -s $si -o $outf.n10.txt

echo "--- $(date) --- Done ---"
