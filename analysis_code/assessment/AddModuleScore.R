
###--------addmodulescore

  #Calculating addModuleScore for each regulon
  DefaultAssay(single_cell) <- "RNA"
  single_cell <- AddModuleScore(single_cell,
                                features = list(Module_regulon),
                                name="Module_regulon")

  add_cell <- as.data.frame(Idents(single_cell))

  add_cell$add_score <- single_cell$Module_regulon1
  colnames(add_cell) <- c("celltypes","regulon_score")
  results_add <- aggregate(regulon_score~celltypes,add_cell,mean)

  results_add$regluons <- Module_regulon[1]

  head(results_add)
