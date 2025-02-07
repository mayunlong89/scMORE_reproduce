##------Figure 7E

Module1 <- regulons[regulons$TF == "DMBX1",]
Module2 <- regulons[regulons$TF == "BCL11A",]
Module3 <- regulons[regulons$TF == "STAT4",]

Module_regulon1 <- c(unique(Module1$TF), unique(Module1$Target))
Module_regulon2 <- c(unique(Module2$TF), unique(Module2$Target))
Module_regulon3 <- c(unique(Module3$TF), unique(Module3$Target))

write.csv(Module_regulon1, file="DMBX1_regulon_target_genes.csv",quote=F, row.names = F)
write.csv(Module_regulon2, file="BCL11A_regulon_target_genes.csv",quote=F, row.names = F)
write.csv(Module_regulon3, file="STAT4_regulon_target_genes.csv",quote=F, row.names = F)

