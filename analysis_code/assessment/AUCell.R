library(AUCell)

cells_rankings <- AUCell_buildRankings(single_cell@assays$RNA@data)

#Calculating AUCell for each regulon
#cells_rankings <- AUCell_buildRankings(organoid_brain_multi@assays$RNA@data)
cells_AUC <- AUCell_calcAUC(Module_regulon,cells_rankings,aucMaxRank = nrow(cells_rankings)*0.05)
#cells_assignment <- AUCell_exploreThresholds(cells_AUC,plotHist = TRUE, nCores = 1, assign =TRUE)

#Extracting AUCell score
AUCell_auc <- as.numeric(getAUC(cells_AUC)[1,])

#Assign AUCell score to each cell
auc_cell <- as.data.frame(Idents(single_cell))

auc_cell$AUCell_auc <- AUCell_auc

colnames(auc_cell) <- c("celltypes","scores")
results_auc <- aggregate(scores~celltypes,auc_cell,mean)

results_auc$regluons <- Module_regulon[1]

head(results_auc)

