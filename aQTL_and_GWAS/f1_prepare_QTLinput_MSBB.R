pheno_path <- "./asapa_relative_usage.rename.FP.filter.txt"
meta_path <- "./specimenID_syn_region_WGS_raceW.txt"    # metadata with rnaseqid to brain region
genotype_PC_path <- "./all_chr.final_rename_FP.PCs.eigenvec" 
stats_path <- "./seq_stat.mat"
output_name <- "asapa_relative_usage.rename.FP"
region <- "FP"

#### prepare phenotype bed/matrix 
args = commandArgs(TRUE)
pheno_path <- args[1]
meta_path <- args[2]
genotype_PC_path <- args[3]
output_name <- args[4] 
stats_path <- args[5]
region <- args[6]

#### 1.1 select RNA samples with corresponding WGS data 
pheno <- read.table(pheno_path, row.names=1, header=T, check.names = F ) 
pheno <- as.matrix(pheno) 
dim(pheno)

# logit transfer and scale 
epsilon <- 1e-5
usage_adj <- pmax(pmin(pheno, 1 - epsilon), epsilon)
pheno <- log(usage_adj / (1 - usage_adj))
pheno <- pheno[apply(pheno, 1, function(x) sd(x, na.rm = TRUE) > 0), ]
pheno_norm <- t(scale(t(pheno))) 

cat("writing:", paste0("./asapa_relative_usage.rename.",region,".filter.norm_scale.txt"),"\n") 
write.table(pheno_norm, paste0("./asapa_relative_usage.rename.",region,".filter.norm_scale.txt"), row.names=T, col.names=T,sep="\t", quote=F) 

#### 1.3 calculate phenotype PCs  
library(PCAForQTL) 

### PC for top 20% various events  
sds <- apply(pheno_norm, 1, sd, na.rm = TRUE)
cutoff <- quantile(sds, 0.8, na.rm = TRUE) 
pheno_norm_top <- pheno_norm[sds >= cutoff, ]
prcompResult <- prcomp(t(pheno_norm_top), center=TRUE, scale.=TRUE) 

PCs <- prcompResult$x 
resultRunElbow <- runElbow(prcompResult = prcompResult) 
print(resultRunElbow) 

pheno_PCs <- as.matrix(PCs[,1:resultRunElbow])
write.table(pheno_PCs, paste0("./asapa_relative_usage.rename.",region,".normalize.wgs_samples.PCs"), quote=F, row.names=T, col.names=T, sep="\t")

#### 1.4 generate tensorQTL input files 
#### phenotype bed 
res <- strsplit(rownames(pheno_norm), split = "[|:-]") 
coord_mat <- do.call(rbind, lapply(res, function(x) x[2:4]))
colnames(coord_mat) <- c("chr", "start", "end")
coord_mat <- data.frame(coord_mat)
coord_mat$event <- rownames(pheno_norm)

pheno_input <- cbind(coord_mat, data.frame(pheno_norm, check.names = F))
write.table(pheno_input, paste0("./asapa_relative_usage.rename.",region,".normalize.wgs_samples.bed"), quote=F, row.names=F, col.names=T, sep="\t")

#### covariants matrix 
library(dplyr)

meta <- read.table(meta_path, header=T, check.names=F, sep="\t")
meta <- meta[meta$tissue == region,]

meta <- meta[, c("specimenID.x", "Braak", "RIN", "sex", "ageDeath", "pmi", "mapped", "rRNA.rate")] 
meta$mapped <- meta$mapped * 1e-6
meta$ageDeath <- gsub("\\+","",meta$ageDeath)

row.names(meta) <- meta$specimenID.x
meta <- meta[, -1]
meta$sex <- as.integer(factor(meta$sex))
meta$pmi <- as.numeric(factor(meta$pmi))
meta$ageDeath <- as.numeric(factor(meta$ageDeath))

genotype_PC <- read.table(genotype_PC_path, row.names=1)
genotype_PC <- genotype_PC[,2:4]
names(genotype_PC) <- c("geno_PC1", "geno_PC2", "geno_PC3")

genotype_PC$sample <- row.names(genotype_PC)
pheno_PCs <- data.frame(pheno_PCs, check.names = F)
pheno_PCs$sample <- row.names(pheno_PCs)
meta$sample <- row.names(meta)

cov_input <- inner_join(meta, genotype_PC, by = "sample")
cov_input <- inner_join(cov_input, pheno_PCs, by = "sample")
rownames(cov_input) <- cov_input$sample
cov_input <- cov_input %>% select(-sample)
cov_input <- t(cov_input)

cov_input <- cov_input[, names(pheno_input)[-(1:4)]]

write.table(cov_input, paste0("./asapa_relative_usage.rename.",region,".normalize.wgs_samples.braak.covariants"), quote=F, row.names=T, col.names=T, sep="\t")







