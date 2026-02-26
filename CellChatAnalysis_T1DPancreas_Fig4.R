
## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/CellChat")
remove.packages("Seurat")
install.packages("https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.3.0.tar.gz", repos = NULL, type = "source")
packageVersion("Seurat")
remove.packages("SeuratObject")
install.packages("https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.1.3.tar.gz", repos = NULL, type = "source")

library(Seurat)
library(SoupX)
library(parallel)
library(purrr)
library(tibble)
library(presto)
library(dplyr)
library(patchwork)
library(plyr)
library(RColorBrewer)
library(multtest)
library(metap) 
library(ggprism)
BiocManager::install("BiocNeighbors")
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
library(glmGamPoi)
devtools::install_github("jinworks/CellChat")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
devtools::install_github('immunogenomics/presto')
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T) 

## ----load Object and subset objects-----
NOD_T1D_Timepoints<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v8.rds")
unique(NOD_T1D_Timepoints$group)
table(NOD_T1D_Timepoints$sample,NOD_T1D_Timepoints$group)


DimPlot(NOD_T1D_Timepoints)
Idents(NOD_T1D_Timepoints)<-NOD_T1D_Timepoints$CellSubType
NOD_T1D_Timepoints$samples<-NOD_T1D_Timepoints$sample
NOD_T1D_Timepoints$sample<-NULL
Progressor <- subset(NOD_T1D_Timepoints, subset = group == "Progressor")
NonProgressor <- subset(NOD_T1D_Timepoints, subset =  group == "Non-Progressor")

# Single Datasets ####
## ----Progressor-----


cellChat_Progressor  <- createCellChat(object = Progressor, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- (CellChatDB)

# set the used database in the object
cellChat_Progressor@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellChat_Progressor <- subsetData(cellChat_Progressor) # This step is necessary even if using the whole database

library(future)
options(future.globals.maxSize = 5e9)  # Set limit to 1GB (adjust based on available memory)
future::plan("multisession", workers = 7) # do parallel
cellChat_Progressor <- identifyOverExpressedGenes(cellChat_Progressor)
cellChat_Progressor <- identifyOverExpressedInteractions(cellChat_Progressor)

gc()  # Free up memory before running the next computation

ptm = Sys.time()
cellChat_Progressor <- computeCommunProb(cellChat_Progressor, type = "triMean")
cellChat_Progressor <- filterCommunication(cellChat_Progressor, min.cells = 10)
cellChat_Progressor <- computeCommunProbPathway(cellChat_Progressor)
cellChat_Progressor <- aggregateNet(cellChat_Progressor)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
ptm = Sys.time()
groupSize <- as.numeric(table(cellChat_Progressor@idents))

par(mfrow = c(1,1), xpd=TRUE)
#netVisual_circle(cellChat_W6_Progressor@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat_Progressor@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength-Week 6 Progressor")

# Compute the network centrality scores
cellChat_Progressor <- netAnalysis_computeCentrality(cellChat_Progressor, slot.name = "netP") 

getwd()
saveRDS(cellChat_Progressor, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/CellChat/cellChat_Progressor.rds")
sessionInfo()

## ----NonProgressor-----
unique(NonProgressor$samples)
### 1. Create cellchat object
cellChat_NonProgressor  <- createCellChat(object = NonProgressor, group.by = "ident", assay = "RNA")

#Set Ligand Receptor Database
CellChatDB <- CellChatDB.mouse 
CellChatDB.use <- (CellChatDB)

# set the used database in the object
cellChat_NonProgressor@DB <- CellChatDB.use
options(future.globals.maxSize = 1e9)  # Set limit to 1GB (adjust based on available memory)
future::plan("multisession", workers = 7) 
### 2. Preprocessing the expression data
# subset the expression data of signaling genes for saving computation cost
cellChat_NonProgressor <- subsetData(cellChat_NonProgressor) # This step is necessary even if using the whole database

cellChat_NonProgressor <- identifyOverExpressedGenes(cellChat_NonProgressor)
cellChat_NonProgressor <- identifyOverExpressedInteractions(cellChat_NonProgressor)

gc()  # Free up memory before running the next computation

### 3. Compute the communication probability
ptm = Sys.time()
cellChat_NonProgressor <- computeCommunProb(cellChat_NonProgressor, type = "triMean")
cellChat_NonProgressor <- filterCommunication(cellChat_NonProgressor, min.cells = 10)

### 4. Infer communication at signaling pathway level

cellChat_NonProgressor <- computeCommunProbPathway(cellChat_NonProgressor)

### 5. Calculate the aggregated cell-cell communication network
cellChat_NonProgressor <- aggregateNet(cellChat_NonProgressor)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(cellChat_NonProgressor@idents))
par(mfrow = c(1,1), xpd=TRUE)
#netVisual_circle(cellChat_W6_NonProgressor@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat_NonProgressor@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength- Non-Progressor")

# Compute the network centrality scores
cellChat_NonProgressor <- netAnalysis_computeCentrality(cellChat_NonProgressor, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

saveRDS(cellChat_NonProgressor, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/CellChat/cellChat_NonProgressor.rds")
sessionInfo()


# Multiple Datasets Comparison ####

## Progressor vs Non-Progressor ----
getwd()
cellChat_NonProgressor <- readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/CellChat/cellChat_NonProgressor.rds")
cellChat_NonProgressor<-updateCellChat(W6_NonProgressor)
cellChat_Progressor <- readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/CellChat/cellChat_Progressor.rds")
cellChat_Progressor<-updateCellChat(cellChat_Progressor)
unique(cellChat_NonProgressor@idents)
unique(cellChat_Progressor@idents)
# Define the cell labels to lift up by combining both cell labels from the conditions

object.list <- list(NonProgressor = cellChat_NonProgressor, Progressor = cellChat_Progressor)
cellchat_PvsNP <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

### 1. Identify altered interactions and cell populations 


par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat_PvsNP, weight.scale = T, measure = "weight", vertex.size.max = 5,vertex.label.cex = 2.2,top = 0.25)

gg1 <- netVisual_heatmap(cellchat_PvsNP, measure = "weight",font.size = 12, font.size.title = 18)
#> Do heatmap based on a merged object
gg1

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg3 <- list()
for (i in 1:length(object.list)) {
  gg3[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+ xlim(0,8) + ylim(0,9)
}
patchwork::wrap_plots(plots = gg3)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_PvsNP, idents.use = "Macrophage")
gg1 + ggrepel::geom_text_repel(aes(label = labels), size = 1.1)

unique(cellchat_PvsNP@meta$ident)
### 2. Identify altered signaling with distinct interaction strength 
gg1 <- rankNet(cellchat_PvsNP, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_PvsNP, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

### 3.Compare outgoing (or incoming) signaling patterns associated with each cell population 
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
# pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 18,font.size = 11,font.size.title = 14)
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 18,font.size = 11,font.size.title = 14)
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# 
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 18, color.heatmap = "GnBu",font.size = 11,font.size.title = 14)
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 18, color.heatmap = "GnBu",font.size = 11,font.size.title = 14)
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 10, height = 18, color.heatmap = "OrRd",font.size = 11,font.size.title = 14)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 10, height = 18, color.heatmap = "OrRd",font.size = 11,font.size.title = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

### 4.Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat_PvsNP, sources.use = 2, targets.use = c(1,3:17),  comparison = c(1, 2), angle.x = 45)


saveRDS(cellchat_PvsNP, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/CellChat/cellchat_PvsNP.rds")

Blue  : "#3A6EA5"
  Green : "#3C8D5A"
  Red   : "#B2182B"
