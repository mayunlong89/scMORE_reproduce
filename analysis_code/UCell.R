#---------UCell
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("UCell")


  library(UCell)

  seu_matrix <- single_cell@assays$RNA@data
  gene.sets <- list(Module_regulon)

  UCellscores <- ScoreSignatures_UCell(seu_matrix, features=gene.sets)
  head(UCellscores)

  results_ucell <- as.data.frame(Idents(single_cell))

  results_ucell$ucell_score <-  UCellscores
  colnames(results_ucell) <- c("celltypes","regulon_score")
  results_ucell2 <- aggregate(regulon_score~celltypes,results_ucell,mean)

  results_ucell2$regluons <- Module_regulon[1]
  colnames(results_ucell2) <- c("celltypes","scores","regulons")
  head(results_ucell2)
