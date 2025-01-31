##-----------Figure 5G

setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/2024-11-20-ctDRTF_re-analysis/10-8psychiatric_disorders/")

grn_outputs <- readRDS("../10_organoid_brain_scHOB/Brain_pando_grn_outputs.rds")

regulons <- grn_outputs$grn
head(regulons)

Module <- regulons[regulons$TF == "GLI3",]
Module_regulon <- c(unique(Module$TF), unique(Module$Target))

write.csv(Module_regulon, file="GLI3_regulon_target_genes.csv",quote=F, row.names = F)

