
mat <- read.table("./asapa_relative_usage.txt", sep = "\t", header = T, row.names = 1, check.names = F)
meta <- read.table("./specimenID_syn_region_WGS_raceW.txt", sep = "\t") # metadata with rnaseqid to brain region

for (region in c("FP", "IFG", "STG", "PHG")){
    ##### extract matrix for each brain region 
	meta2 <- meta[meta$V3 == region,]
    
    old_names <- colnames(mat)
    id <- sub("\\.rmdup$", "", old_names)
    idx <- match(id, meta2$V2)
    keep <- !is.na(idx)
    mat2 <- mat[, keep, drop = FALSE]
    colnames(mat2) <- meta2$V1[idx[keep]]

	print(dim(mat2))
    cat("writing:", paste0("./asapa_relative_usage.rename.",region,".txt"),"\n") 
    write.table(mat2, paste0("./asapa_relative_usage.rename.",region,".txt"), quote = F, row.names = T, col.names = T, sep = "\t") 

    #### filter, logit-transfer, and scale matrix 
    ##### filter matrix 
    impute_row_mean <- function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
    }

    # 0. select exon node with usage > 5% in at least 10% samples 
    data_filtered1 <- mat2[rowSums(mat2 >= 0.05, na.rm = TRUE) >= 0.1*ncol(mat2),]
    cat("filter_step1_remain: ", nrow(data_filtered1),"\n") 
    # 1. remove NA > 90% 
    na_threshold <- 0.9 * ncol(data_filtered1)
    data_filtered2 <- data_filtered1[rowSums(is.na(data_filtered1)) <= na_threshold, ] 
    cat("filter_step2_remain: ", nrow(data_filtered2),"\n") 
    # 3. impute mean value 
    data_imputed <- t(apply(data_filtered2, 1, impute_row_mean)) 

    # 4. sd > 0 
    data_imputed_scaled <- t(scale(t(data_imputed))) 
    sd <- apply(data_imputed_scaled, 1, sd) 
    keep <- sd > 0  
    data_filtered3 <- na.omit(data_imputed[keep, ]) 
    cat("filter_step3_remain: ", nrow(data_filtered3),"\n")

	print(dim(data_filtered3))
    cat("writing:", paste0("./asapa_relative_usage.rename.",region,".filter.txt"),"\n") 
    write.table(data_filtered3, paste0("./asapa_relative_usage.rename.",region,".filter.txt"), row.names=T, col.names=T,sep="\t", quote=F) 
}




