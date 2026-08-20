#! /usr/bin/bash
#$ -cwd
#$ -N apaQTL
#$ -o loga/job.$JOB_NAME.$TASK_ID.out
#$ -e loga/job.$JOB_NAME.$TASK_ID.err
#$ -V
#$ -l h_rt=8:00:00,h_data=8G
#$ -pe shared 1
#$ -t 1-1:1

idx=$SGE_TASK_ID

source $HOME/.bash_profile

set -x -e -o pipefail
. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO

####################

DIR_BASE=/home/directory
DIR_WORK=${DIR_BASE}/asAPA_analysis
DIR_RESULTS=${DIR_WORK}/aQTL/results
DIR_PLOTS=${DIR_RESULTS}/plots
DIR_PLOTS_CONC=${DIR_RESULTS}/plots_conc
mkdir -p ${DIR_PLOTS}

FP=${DIR_RESULTS}/asapa_relative_usage.rename.FP.top.cis_qtl.txt
IFG=${DIR_RESULTS}/asapa_relative_usage.rename.IFG.top.cis_qtl.txt
PG=${DIR_RESULTS}/asapa_relative_usage.rename.PHG.top.cis_qtl.txt
STG=${DIR_RESULTS}/asapa_relative_usage.rename.STG.top.cis_qtl.txt

# List of asAPA SNPs and genes
# SNPs: chrom | pos | ref | alt
# Genes: ENSG ID | Gene name
asapa_snps=${DIR_RESULTS}/asapa_snps.txt
control_snps=${DIR_RESULTS}/control_snps.txt
asapa_genes=${DIR_RESULTS}/asapa_genes.txt
control_genes=${DIR_RESULTS}/control_genes.txt

# SNPs with allelic ratios for reference allele
# SNPs: chrom | pos | ref | alt | ref count | alt count | ref ratio
asapa_ratio=${DIR_RESULTS}/asapa_snps_ratio_ref.txt

conc_snps=${DIR_RESULTS}/conc_snps.txt
conc_genes=${DIR_RESULTS}/conc_genes.txt

echo "--- $(date) --- Make Plots ---"

Rscript ${DIR_WORK}/scripts_review/apa_QTL_plots.R \
    ${FP} ${IFG} ${PG} ${STG} \
    ${asapa_snps} \
    ${control_snps} \
    ${asapa_genes} \
    ${control_genes} \
    ${DIR_PLOTS} \
    ${asapa_ratio}

Rscript ${DIR_WORK}/scripts_review/apa_QTL_plots.R \
    ${FP} ${IFG} ${PG} ${STG} \
    ${conc_snps} \
    ${control_snps} \
    ${conc_genes} \
    ${control_genes} \
    ${DIR_PLOTS_CONC} \
    ${asapa_ratio}

rm -f ./VennDiagram.*.log

echo "--- $(date) --- Done ---"
