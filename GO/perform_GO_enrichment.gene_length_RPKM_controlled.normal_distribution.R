#converts inputted lists of genes from gene symbol annotation to Ensembl IDs
#performs GO enrichment using a permutation test, controlling for gene length and gene expression (tpm)
#uses normal ditribution rather than nonparametric distribution to test
#(adapted from Stephen's script)
rm(list=ls())

#GLOBAL VARIABLES
#change as desired

#GENE_LENGTH_MATCH_PERCENT = 0.05  #gene length of matched controls must be within this fraction of query gene length
GENE_LENGTH_MATCH_PERCENT = .20 #0.10  #gene length of matched controls must be within this fraction of query gene length
NUMBER_OF_PERMUTATIONS = 10000  #number of permutations in GO tests
#GENE_tpm_MATCH_PERCENT = 0.1 #gene tpm value of matched controls must be within this fraction of query gene tpm
GENE_TPM_MATCH_PERCENT = .20 #0.10 #gene tpm value of matched controls must be within this fraction of query gene tpm
#THROW_OUT_GENES = "yes" # yes or no, whether or not to throw out query genes missing unmatched controls 
########################################

args = commandArgs(trailingOnly=TRUE)

tpm_file = '/path/to/TPM/average_tpm.txt'
#tpm_file = args[2]
gene_annotation_fn = '/path/to/gene_annotation.bed'

#args = commandArgs(trailingOnly=TRUE)
#QUERY_STATES = args[1]
#THROW_OUT_GENES = args[2]

library(biomaRt)
library(data.table)
library(dplyr)
#library(data.table)
#biomartCacheClear()
unlink(tools::R_user_dir("biomaRt", which = "cache"), recursive = TRUE)
#ensembl = useMart( biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
#load('../results/gene_label_after_GENCODE_CNV_filtering.RData')
#all_asrs_anno = fread('../snv_types/snv_gene_anno_AR/all_cell_lines.asRS.gene.bed')
#input_query_list_fn = 'GO_queries/Systemic_lupus_erythematosus_ensembl_genes.txt'
#input_query_list_fn = args[1] 
#background_list_fn = args[2]
#background_GO_table_fn = args[3]
#output_file = args[4]
#input_query_list_fn = '../results/all_EI_AEI_ensembl_ETG_sans_Jae.txt'

input_query_list_fn = args[1]
prefix = gsub(".txt", "", basename(input_query_list_fn))
prefix = paste0(prefix, '_gene_length_match_', GENE_LENGTH_MATCH_PERCENT, '_gene_TPM_match_', GENE_TPM_MATCH_PERCENT)
GO_query_table_fn = paste0('../GO/', prefix, '_query_GO_table_gene_length_match.RData')
background_GO_table_fn = paste0('../GO/', prefix, '_background_GO_table.RData')
output_file = paste0('../GO/', prefix, 'GO_enrichment_res.txt')
#GO_query_table_fn = paste0('../results/all_universal_ensembl_ETG_query_GO_table_gene_length_RPKM_controlled.RData')
#background_GO_table_fn = '../results/all_universal_ensembl_ETG_background_GO_table_gene_length_RPKM_controlled.RData'
#output_file = '../results/all_universal_ensembl_ETG.GO_enrichment_res_gene_length_RPKM_controlled.txt'





GENE_LENGTH_MATCH_PERCENT = as.numeric(GENE_LENGTH_MATCH_PERCENT)
GENE_TPM_MATCH_PERCENT = as.numeric(GENE_TPM_MATCH_PERCENT)

input_query_list = fread(input_query_list_fn)
#gene_annotation = fread(gene_annotation_fn)




tpm_df = read.table(tpm_file, header=TRUE, sep='\t',stringsAsFactors=FALSE)
#tpm_df = fread(tpm_file)
colnames(tpm_df) = c("Gene", "tpm")
#background_list = scan(background_list,what='')
#the background list is the tpm table gene names
background_list = filter(input_query_list, label == 'testable')
#input_query_list = scan(input_query_list,what='')
input_query_list = filter(input_query_list, label != 'testable')
input_query_list =  unique(input_query_list$ensembl_id)
background_list = setdiff(background_list$ensembl_id, input_query_list)


library(dplyr)
library(data.table)
library(biomaRt)
#biomartCacheClear()
unlink(tools::R_user_dir("biomaRt", which = "cache"), recursive = TRUE)
#ensembl = useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')

ensembl=useMart("ensembl", dataset='hsapiens_gene_ensembl')

#setEnsemblSSL(list(ssl_verifypeer = FALSE))
#httr::set_config(httr::config(ssl_verifypeer = FALSE))
#ensembl = useEnsembl( biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
cat('getting matched gene length background list\n')
gene_annotation = fread(gene_annotation_fn)
gene_annotation$length = gene_annotation$V3 - gene_annotation$V2
gene_annotation = gene_annotation %>% rename(ensembl_id = V4)

query_gene_lengths = filter(gene_annotation, ensembl_id %in% input_query_list)

background_gene_lengths = filter(gene_annotation, ensembl_id %in% background_list)




# some query genes are not in tpm_df, prob low coverage
table(query_gene_lengths$ensembl_gene_id %in% tpm_df$Gene)
query_gene_lengths = filter(query_gene_lengths, ensembl_id %in% tpm_df$Gene)


#now let's get matched gene lists
matched_gene_list = list()
number_genes_match_length = rep(0, nrow(query_gene_lengths))
number_genes_match_tpm = rep(0, nrow(query_gene_lengths))
for(i in 1:nrow(query_gene_lengths)){
	if(i %% 100 == 0) print(i)
	query_gene_name = query_gene_lengths$ensembl_id[i]
	query_gene_length = query_gene_lengths$length[i]
	#if (query_gene_name %in% background_list) {
	  query_gene_tpm = tpm_df[which(tpm_df$Gene == query_gene_name),]$tpm 
	  lower_bound = query_gene_length - query_gene_length* GENE_LENGTH_MATCH_PERCENT
	  upper_bound = query_gene_length + query_gene_length* GENE_LENGTH_MATCH_PERCENT
	  tpm_lower_bound = query_gene_tpm - query_gene_tpm* GENE_TPM_MATCH_PERCENT
	  tpm_upper_bound = query_gene_tpm + query_gene_tpm* GENE_TPM_MATCH_PERCENT

	  if(tpm_upper_bound < tpm_lower_bound){
	  print('warning: tpm upper is lower than tpm lower')
	  }

	  background_rows_to_keep_bylength = background_gene_lengths$length >= lower_bound & background_gene_lengths$length <= upper_bound
	  background_genes_bylength = background_gene_lengths[background_rows_to_keep_bylength,'ensembl_id']
	  number_genes_match_length[i] = sum(background_rows_to_keep_bylength)

	  background_genes_bylength_tpm = tpm_df[tpm_df$Gene %in% background_genes_bylength$ensembl_id,]
	  background_rows_to_keep_bylengthxtpm = background_genes_bylength_tpm$tpm >= tpm_lower_bound & background_genes_bylength_tpm$tpm <= tpm_upper_bound
	  background_genes = background_genes_bylength_tpm[background_rows_to_keep_bylengthxtpm,'Gene']

	  matched_gene_list[[query_gene_name]] = background_genes
	#}
}

#now make sure we background genes for each query gene
number_background_matched_list = sapply(matched_gene_list,function(i) length(i))
table(number_background_matched_list)
#if(any(number_background_matched_list == 0)){
#	cat('some of the query genes do not have any matched gene length controls. Please increase you GENE_LENGTH_MATCH_PERCENT global variable value\n')
#	quit()
#} #query_gene_lengths[number_background_matched_list ==0,]
#cat('count number of query genes without any matched controls\n')
throw_out_queries = query_gene_lengths$ensembl_id[number_background_matched_list ==0]
#print(filter(gene_annotation, ensembl_id %in% throw_out_queries))
for(bad_gene in throw_out_queries){
	matched_gene_list[[bad_gene]] = NULL
}
print(paste(throw_out_queries,collapse=","))

## ELAINE UNCOMMENTED THIS LINE 12/06/2023 BECAUSE THE NUMBER OF QUERY GENES WAS TOO DIFFERENT FROM NUMBER OF BACKGROUND GENES 
#if(THROW_OUT_GENES == "yes"){
##query_gene_lengths = query_gene_lengths[-c(which(query_gene_lengths$ensembl_gene_id %in% throw_out_queries)),] ##Mudra added this line to throw out the query genes that do not have matched background genes. However, sometimes this does not yield good GO results.
#cat(paste('number of query genes thrown out: ',length(throw_out_queries),'\n',sep=""))
#} else if (THROW_OUT_GENES == "no") {
#cat("keep unmatched query genes")
#} else {
#cat("whether to throw out unmatched query genes is not specified, keep them by default")
#}


total_genes = length(query_gene_lengths$ensembl_id)

#next get GO terms
cat('getting GO terms\n')
packageVersion("dbplyr")
packageVersion("dplyr")
unlink(tools::R_user_dir("biomaRt", which = "cache"), recursive = TRUE)
background_GO_table = getBM(attributes=c('hgnc_symbol','ensembl_gene_id','go_id','name_1006','definition_1006','namespace_1003'),filters=c("ensembl_gene_id"),values=substr(background_gene_lengths$ensembl_id,1,15),mart=ensembl)
query_GO_table = getBM(attributes=c('hgnc_symbol','ensembl_gene_id','go_id','name_1006','definition_1006','namespace_1003'),filters=c("ensembl_gene_id"),values=substr(query_gene_lengths$ensembl_id,1,15), mart=ensembl) ##change this so that you remove the "bad genes" from the query gene's 

background_GO_table = filter(background_GO_table, !ensembl_gene_id %in% substr(throw_out_queries, 1, 15))
save(query_GO_table, file = GO_query_table_fn)
save(background_GO_table, file = background_GO_table_fn)
#save(matched_gene_list, file = matched_gene_list_fn)


#remove '' go_id
query_GO_table = query_GO_table[query_GO_table$go_id != '',]
query_GO_terms = query_GO_table$go_id #these are the query GO terms
unique_query_GO_terms = query_GO_terms[!duplicated(query_GO_terms)]
unique_query_GO_table = query_GO_table[!duplicated(query_GO_terms),]



#cat('saving intermediate .RData session in case permutation test is interrupted\n')
cat('running permutations for GO enrichment\n')
#run permutations 10,000 times. makes sure to store the GO terms for the background list for teach of the 10,000 times. This is important for the sake of not having to run this test 10,000 times for each GO term.
#let's try doing this in parallel for speed up
library(doParallel)
library(foreach)
#noCores = detectCores() -1
noCores = 5 #detectCores() was giving a very large number
cl = makeCluster(noCores,outfile="") #outfile ='' to allow print statements to write to master node
registerDoParallel(cl,cores=noCores)

GO_terms_per_permutation = list()
#for(i in 1:NUMBER_OF_PERMUTATIONS){
system.time({
GO_terms_per_permutation = foreach(i = 1:NUMBER_OF_PERMUTATIONS) %dopar%{
	if (i %% 1000 ==0)print(i)
	#get random sample
	random_samples = sapply(matched_gene_list,function(j)sample(j,1))		
	random_samples = substr(random_samples, 1, 15)
#debugging
	GO_ids = background_GO_table[background_GO_table$ensembl_gene_id %in% random_samples,'go_id']
	#remove the "" GO_IDs (which should be here
	GO_ids = GO_ids[GO_ids != '']
	#remove go_ids that are not in the query dataset. These are never going to be tested
	GO_ids = GO_ids[GO_ids %in% unique_query_GO_terms]
	return(GO_ids)
	#GO_terms_per_permutation[[i]] =  GO_ids #store the go_id terms 
}
})



#next get enrichment for each GO term represented in the query genes
cat('getting GO enrichment\n')
#okay now let's get enrichment for each GO term
#make sure you're using max number of processors available using getDoParWorkers() function
#for(i in 1:length(unique_query_GO_terms)){
#this next part takes a long time to start for some reason. But once it starts it goes pretty quickly
system.time({
tmp = foreach(i=1:length(unique_query_GO_terms),.combine='rbind') %dopar%{
#tmp = foreach(i=1:20,.combine='rbind') %dopar%{
	if(i %% 50 ==0)print(i)
	#print(i)
	the_GO_term = unique_query_GO_terms[i]
	#print(the_GO_term)
	the_ocurrences_in_query = sum(query_GO_terms == the_GO_term)
	#ocurrences_in_query[i] = sum(query_GO_terms == the_GO_term)
	ocurrences_in_background = sapply(GO_terms_per_permutation,function(j)sum(j==the_GO_term))
	#enrichment_p_values[i] = sum(ocurrences_in_background >= ocurrences_in_query[i])/NUMBER_OF_PERMUTATIONS
	#the_enrichment_p_value = sum(ocurrences_in_background >= the_ocurrences_in_query)/NUMBER_OF_PERMUTATIONS
	the_enrichment_p_value = pnorm(the_ocurrences_in_query, mean = mean(ocurrences_in_background), sd= sd(ocurrences_in_background), lower.tail = FALSE, log.p=F)
#	the_enrichment_p_value = 1
	return(c(the_ocurrences_in_query,the_enrichment_p_value))
}
})

ocurrences_in_query = tmp[,1]
enrichment_p_values = tmp[,2]
FDR = p.adjust(enrichment_p_values,'fdr')

stopCluster(cl) #to stop using the cluster you created

output_table = data.frame(unique_query_GO_terms,ocurrences_in_query,enrichment_p_values,FDR)
output_table = merge(output_table,unique_query_GO_table,by.x='unique_query_GO_terms',by.y='go_id',all.x=TRUE,all.y=FALSE)
output_table$hgnc_symbol = NULL
output_table$ensembl_gene_id = NULL

#reorder table by p values
output_table = output_table[order(output_table$enrichment_p_values),]

output_table$total_genes = total_genes

write.table(output_table,file=output_file,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

#cat('saving results\n')
#save.image()

cat('job completed\n')
