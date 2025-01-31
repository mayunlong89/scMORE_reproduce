
##---------Supplemental--Figure S24
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/10-8psychiatric_disorders/")

grn_outputs <- readRDS("../10_organoid_brain_scHOB/Brain_pando_grn_outputs.rds")

regulons <- grn_outputs$grn
head(regulons)



Module <- regulons[regulons$TF == "PLAGL1",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="PLAGL1_regulon_target_genes.csv",quote=F, row.names = F)




Module <- regulons[regulons$TF == "PRDM12",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="PRDM12_regulon_target_genes.csv",quote=F, row.names = F)



Module <- regulons[regulons$TF == "SOX3",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="SOX3_regulon_target_genes.csv",quote=F, row.names = F)


Module <- regulons[regulons$TF == "PRDM13",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="PRDM13_regulon_target_genes.csv",quote=F, row.names = F)




