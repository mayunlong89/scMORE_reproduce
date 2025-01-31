
##---------Supplemental--Figure S23
setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/10-8psychiatric_disorders/")

grn_outputs <- readRDS("../10_organoid_brain_scHOB/Brain_pando_grn_outputs.rds")

regulons <- grn_outputs$grn
head(regulons)

##ExN-specific neurons----------------------------------------------------------
Module <- regulons[regulons$TF == "NEUROD6",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="NEUROD6_regulon_target_genes.csv",quote=F, row.names = F)


Module <- regulons[regulons$TF == "MEF2C",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="MEF2C_regulon_target_genes.csv",quote=F, row.names = F)


Module <- regulons[regulons$TF == "BHLHE22",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="BHLHE22_regulon_target_genes.csv",quote=F, row.names = F)


Module <- regulons[regulons$TF == "POU2F2",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="POU2F2_regulon_target_genes.csv",quote=F, row.names = F)



Module <- regulons[regulons$TF == "TFAP2C",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="TFAP2C_regulon_target_genes.csv",quote=F, row.names = F)



Module <- regulons[regulons$TF == "TBR1",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="TBR1_regulon_target_genes.csv",quote=F, row.names = F)



Module <- regulons[regulons$TF == "NPAS4",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="NPAS4_regulon_target_genes.csv",quote=F, row.names = F)

Module <- regulons[regulons$TF == "KLF9",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="KLF9_regulon_target_genes.csv",quote=F, row.names = F)

##InN-specific neurons----------------------------------------------------------
Module <- regulons[regulons$TF == "ASCL1",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="ASCL1_regulon_target_genes.csv",quote=F, row.names = F)

Module <- regulons[regulons$TF == "TSHZ2",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="TSHZ2_regulon_target_genes.csv",quote=F, row.names = F)


Module <- regulons[regulons$TF == "MXD3",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="MXD3_regulon_target_genes.csv",quote=F, row.names = F)



Module <- regulons[regulons$TF == "HES4",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="HES4_regulon_target_genes.csv",quote=F, row.names = F)


Module <- regulons[regulons$TF == "ZEB2",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="ZEB2_regulon_target_genes.csv",quote=F, row.names = F)

