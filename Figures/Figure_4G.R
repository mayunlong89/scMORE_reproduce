
##-----------Figure 4G------
####-------------------scMORE monocyte subtypes inference analysis

##-----IBD
gene_info_IBD <- read.table("IBD_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_IBD <- read.csv("IBD_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

gene_info<- gene_info_IBD 
snp_info <- snp_info_IBD


results_IBD_mono <- scMore(single_cell_mono,
                      snp_info,
                      gene_info,
                      n_targets = 5,
                      perm_n = 1000,
                      theta = 0.5,
                      alpha = 1,
                      buffer = 500,
                      top_n = 5,
                      p1 = 0.05,
                      p2 = 0.05,
                      p3 = 0.05,
                      peak2gene_method = 'Signac',
                      infer_method = 'glm',
                      method = 'cosine',
                      nSeed = 1234)


saveRDS(results_IBD_mono, file="results_IBD_scMORE_mono.rds")

#E-statistics
Escore_IBD <- getEnergyScore(results_IBD_mono$scMore_trait_results,targetCelltype = 2)
log2(Escore_IBD+1)+1


write.csv(results_IBD_mono$scMore_trait_results, file="results_IBD_mono_scMore_trait_results.csv",quote=F, row.names = F)


##-----PBC
gene_info_PBC <- read.table("PBC_count_GWAS_final.hg19_SNP_Gene_Analysis_results.genes.out",header = TRUE)
#snp_info
snp_info_PBC <- read.table("PBC_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_PBC
snp_info <- snp_info_PBC


results_PBC_mono <- scMore(single_cell_mono,
                      snp_info,
                      gene_info,
                      n_targets = 5,
                      perm_n = 1000,
                      theta = 0.5,
                      alpha = 1,
                      buffer = 500,
                      top_n = 5,
                      p1 = 0.05,
                      p2 = 0.05,
                      p3 = 0.05,
                      peak2gene_method = 'Signac',
                      infer_method = 'glm',
                      method = 'cosine',
                      nSeed = 1234)


saveRDS(results_PBC_mono, file="results_PBC_scMORE_mono.rds")

#E-statistics
Escore_PBC_mono <- getEnergyScore(results_PBC_mono$scMore_trait_results,targetCelltype = 2)
log2(Escore_PBC_mono+1)+1


write.csv(results_PBC_mono$scMore_trait_results, file="results_PBC_mono_scMore_trait_results.csv",quote=F, row.names = F)



##-----RA
gene_info_RA <- read.table("RA_count_GWAS_final.hg19_SNP_Gene_Analysis_results.genes.out",header = TRUE)
#snp_info
snp_info_RA <- read.table("RA_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_RA
snp_info <- snp_info_RA


results_RA_mono <- scMore(single_cell_mono,
                     snp_info,
                     gene_info,
                     n_targets = 5,
                     perm_n = 1000,
                     theta = 0.5,
                     alpha = 1,
                     buffer = 500,
                     top_n = 5,
                     p1 = 0.05,
                     p2 = 0.05,
                     p3 = 0.05,
                     peak2gene_method = 'Signac',
                     infer_method = 'glm',
                     method = 'cosine',
                     nSeed = 1234)

saveRDS(results_RA_mono, file="results_RA_scMORE_mono.rds")

#E-statistics
Escore_RA <- getEnergyScore(results_RA_mono$scMore_trait_results,targetCelltype = 2)
log2(Escore_RA+1)+1

write.csv(results_RA_mono$scMore_trait_results, file="results_RA_mono_scMore_trait_results.csv",quote=F, row.names = F)



##-----SLE
gene_info_SLE <- read.table("SLE_count_GWAS_final.hg19_SNP_Gene_Analysis_results.genes.out",header = TRUE)
#snp_info
snp_info_SLE <- read.table("SLE_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_SLE
snp_info <- snp_info_SLE


results_SLE_mono <- scMore(single_cell_mono,
                      snp_info,
                      gene_info,
                      n_targets = 5,
                      perm_n = 1000,
                      theta = 0.5,
                      alpha = 1,
                      buffer = 500,
                      top_n = 5,
                      p1 = 0.05,
                      p2 = 0.05,
                      p3 = 0.05,
                      peak2gene_method = 'Signac',
                      infer_method = 'glm',
                      method = 'cosine',
                      nSeed = 1234)


saveRDS(results_SLE_mono, file="results_SLE_scMORE_mono.rds")

#E-statistics
Escore_SLE_mono <- getEnergyScore(results_SLE_mono$scMore_trait_results,targetCelltype = 2)
log2(Escore_SLE_mono+1)+1

write.csv(results_SLE_mono$scMore_trait_results, file="results_SLE_mono_scMore_trait_results.csv",quote=F, row.names = F)


##-----UC
gene_info_UC <- read.table("UC_processed_magma_results.genes.out",header = TRUE)
#snp_info
snp_info_UC <- read.table("UC_maf0.01.txt",header=T,stringsAsFactors = FALSE)

gene_info<- gene_info_UC
snp_info <- snp_info_UC



results_UC_mono <- scMore(single_cell_mono,
                     snp_info,
                     gene_info,
                     n_targets = 5,
                     perm_n = 1000,
                     theta = 0.5,
                     alpha = 1,
                     buffer = 500,
                     top_n = 5,
                     p1 = 0.05,
                     p2 = 0.05,
                     p3 = 0.05,
                     peak2gene_method = 'Signac',
                     infer_method = 'glm',
                     method = 'cosine',
                     nSeed = 1234)

saveRDS(results_UC_mono, file="results_UC_scMORE_mono.rds")

#E-statistics
Escore_UC <- getEnergyScore(results_UC_mono$scMore_trait_results,targetCelltype = 2)
log2(Escore_UC+1)+1

write.csv(results_UC_mono$scMore_trait_results, file="results_UC_mono_scMore_trait_results.csv",quote=F, row.names = F)




##----Heatmap for subset of monocytes analysis

#@ Subset monocyte specific regulons-----heatmap---Figure 4G
mono_data <- read.table("02_mono_subset_regulons_fivediseases.txt",header = TRUE)

mono_data1 <- mono_data [,c(-1)]
rownames(mono_data1) <- mono_data [,1]

#mono_data1_RRS <- mono_data1[,c(1,3,5,7,9,11,13,15,17,19)]
#mono_data1_P <- mono_data1[,c(2,4,6,8,10,12,14,16,18,20)]

mono_data1_RRS <- mono_data1[,c(1,5,9,13,17,3,7,11,15,19)]
mono_data1_P <- mono_data1[,c(2,6,10,14,18,4,8,12,16,20)]

heatmap_data <- t(mono_data1_RRS)
p_vals <- t(mono_data1_P)
pearson_heatmap<-pheatmap(heatmap_data,cellwidth =16, cellheight =9,
                          fontsize=8,
                          border_color="black",
                          scale="column",
                          fontsize_row=8, 
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white",'#336699')))(102),
                          cluster_rows = F,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.05 & heatmap_data > 3, "**", ifelse(p_vals<0.05 & heatmap_data>1,"*","")), nrow(p_vals)))


