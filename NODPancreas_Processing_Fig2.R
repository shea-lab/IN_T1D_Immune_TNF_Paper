#' ---
#' title: "R Notebook"
#' author: "Jyotirmoy Roy adapted from Kate Griffin"
#' Date: "`r Sys.Date()`"
#' output:
#'   pdf_document: default
#'   html_notebook: default
#' editor_options: 
#'   markdown: 
#'     wrap: 72
#' ---
#' 
#' This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
#' execute code within the notebook, the results appear beneath the code.
#' 
#' Try executing this chunk by clicking the *Run* button within the chunk
#' or by placing your cursor inside it and pressing *Cmd+Shift+Enter*.
#' 
## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' Add a new chunk by clicking the *Insert Chunk* button on the toolbar or
#' by pressing *Cmd+Option+I*.
#' 
#' When you save the notebook, an HTML file containing the code and output
#' will be saved alongside it (click the *Preview* button or press
#' *Cmd+Shift+K* to preview the HTML file).
#' 
#' The preview shows you a rendered HTML copy of the contents of the
#' editor. Consequently, unlike *Knit*, *Preview* does not run any R code
#' chunks. Instead, the output of the chunk when it was last run in the
#' editor is displayed.
#' 
#' # Load libraries
#' 
## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq")
library(Seurat)
packageVersion("Seurat")
library(SoupX)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scRNAseq")
library(scRNAseq)

devtools::install_github("ayshwaryas/ddqc_R")
library(ddqcR)

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(parallel)
library(purrr)
library(tibble)
install_github('immunogenomics/presto')
library(presto)
library(dplyr)
library(patchwork)
library(plyr)

library(RColorBrewer)
BiocManager::install("multtest")
library(multtest)
library(metap) 
library(ggprism)
BiocManager::install('glmGamPoi')
library(glmGamPoi)
install.packages("tibble")

#' 
## ----source files-------------------------------------------------------------------------------------------------------------------------------------
source("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/SingleCellRNASeq/remove_soup.R")
source("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/SingleCellRNASeq/find_doublets.R")
source("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/SingleCellRNASeq/run_ddqcR.R")


## ----read in 10x data---------------------------------------------------------------------------------------------------------------------------------
data_dir = "11329_JR_pool/"

sample.names = list.files(path = '11329_JR_pool/')
sample.names

folder.names = paste0(sample.names, "/")

sample_dir_list = folder.names




Week18_1 = remove_soup(Week18_1, data_dir, sample_dir_list[1])
Week18_2 =  remove_soup(Week18_2, data_dir, sample_dir_list[9])
Week18_3 =  remove_soup(Week18_3, data_dir, sample_dir_list[10])
Week18_4 =  remove_soup(NWeek18_4, data_dir, sample_dir_list[11])
Week18_5 =  remove_soup(Week18_5, data_dir, sample_dir_list[12])
Week18_6 =  remove_soup(Week18_6, data_dir, sample_dir_list[13])

Week12_1 = remove_soup(Week12_1, data_dir, sample_dir_list[14])
Week12_2 =  remove_soup(Week12_2, data_dir, sample_dir_list[15])
Week12_3 =  remove_soup(Week12_3, data_dir, sample_dir_list[16])
Week12_4 =  remove_soup(NWeek12_4, data_dir, sample_dir_list[2])

Week6_1 = remove_soup(Week6_1, data_dir, sample_dir_list[3])
Week6_2 =  remove_soup(Week6_2, data_dir, sample_dir_list[4])
Week6_3 =  remove_soup(Week6_3, data_dir, sample_dir_list[5])
Week6_4 =  remove_soup(NWeek6_4, data_dir, sample_dir_list[6])
Week6_5 =  remove_soup(Week6_5, data_dir, sample_dir_list[7])
Week6_6 =  remove_soup(Week6_6, data_dir, sample_dir_list[8])

#save.image(file = "Allergy_Katepiepeline_v1.RData")


#' 
#' # ddqcR
#' 
#' run ddqcR with "filter_ddqcr" function params: seurat object, "sample
#' name", data_dir, sample_dir
#' 
## ----run ddqcR----------------------------------------------------------------------------------------------------------------------------------------

sample_list<- c("Week18_1", "Week18_2", "Week18_3", "Week18_4", "Week18_5","Week18_6",
                "Week12_1","Week12_2","Week12_3","Week12_4",
                "Week6_1","Week6_2","Week6_3","Week6_4","Week6_5","Week6_6")

Week18_1 = filter_ddqcr(Week18_1, sample_list[1], sample_dir_list[1])
Week18_2 =  filter_ddqcr(Week18_2, sample_list[2], sample_dir_list[10])
Week18_3 =  filter_ddqcr(Week18_3, sample_list[3], sample_dir_list[10])
Week18_4 =  filter_ddqcr(Week18_4, sample_list[4], sample_dir_list[11])
Week18_5 =  filter_ddqcr(Week18_5, sample_list[5], sample_dir_list[12])
Week18_6 =  filter_ddqcr(Week18_6, sample_list[6], sample_dir_list[13])

Week12_1 = filter_ddqcr(Week12_1, sample_list[7], sample_dir_list[14])
Week12_2 =  filter_ddqcr(Week12_2, sample_list[8], sample_dir_list[15])
Week12_3 =  filter_ddqcr(Week12_3, sample_list[10], sample_dir_list[16])
Week12_4 =  filter_ddqcr(Week12_4, sample_list[10], sample_dir_list[2])

Week6_1 = filter_ddqcr(Week6_1, sample_list, sample_dir_list[3])
Week6_2 =  filter_ddqcr(Week6_2, sample_list, sample_dir_list[4])
Week6_3 =  filter_ddqcr(Week6_3, sample_list, sample_dir_list[5])
Week6_4 =  filter_ddqcr(Week6_4, sample_list, sample_dir_list[6])
Week6_5 =  filter_ddqcr(Week6_5, sample_list, sample_dir_list[7])
Week6_6 =  filter_ddqcr(Week6_6, sample_list, sample_dir_list[8])



#' 
#' # Run Doublet Finder on Samples
#' 
#' run doublet finder with "find_doublets" function params: seu_object,
#' "sample name", data_dir, sample_dir
#' 
#' 
## ----Run Doublet Finder on Samples--------------------------------------------------------------------------------------------------------------------

Week18_1 = find_doublets(Week18_1, sample_list[1], sample_dir_list[1])
Week18_2 =  find_doublets(Week18_2, sample_list[2], sample_dir_list[10])
Week18_3 =  find_doublets(Week18_3, sample_list[3], sample_dir_list[10])
Week18_4 =  find_doublets(Week18_4, sample_list[4], sample_dir_list[11])
Week18_5 =  find_doublets(Week18_5, sample_list[5], sample_dir_list[12])
Week18_6 =  find_doublets(Week18_6, sample_list[6], sample_dir_list[13])

Week12_1 = find_doublets(Week12_1, sample_list[7], sample_dir_list[14])
Week12_2 =  find_doublets(Week12_2, sample_list[8], sample_dir_list[15])
Week12_3 =  find_doublets(Week12_3, sample_list[10], sample_dir_list[16])
Week12_4 =  find_doublets(Week12_4, sample_list[10], sample_dir_list[2])

Week6_1 = find_doublets(Week6_1, sample_list, sample_dir_list[3])
Week6_2 =  find_doublets(Week6_2, sample_list, sample_dir_list[4])
Week6_3 =  find_doublets(Week6_3, sample_list, sample_dir_list[5])
Week6_4 =  find_doublets(Week6_4, sample_list, sample_dir_list[6])
Week6_5 =  find_doublets(Week6_5, sample_list, sample_dir_list[7])
Week6_6 =  find_doublets(Week6_6, sample_list, sample_dir_list[8])




#' 
## ----save filtered seurat objects---------------------------------------------------------------------------------------------------------------------
saveRDS(Week18_1, file = "./filtered_Week18_1.rds")
saveRDS(Week18_2, file = "./filtered_Week18_2.rds")
saveRDS(Week18_3, file = "./filtered_Week18_3.rds")
saveRDS(Week18_4, file = "./filtered_Week18_4.rds")
saveRDS(Week18_5, file = "./filtered_Week18_5.rds")
saveRDS(Week18_6, file = "./filtered_Week18_6.rds")

saveRDS(Week12_1, file = "./filtered_Week12_1.rds")
saveRDS(Week12_2, file = "./filtered_Week12_2.rds")
saveRDS(Week12_3, file = "./filtered_Week12_3.rds")
saveRDS(Week12_4, file = "./filtered_Week12_4.rds")

saveRDS(Week6_1, file = "./filtered_Week6_1.rds")
saveRDS(Week6_2, file = "./filtered_Week6_2.rds")
saveRDS(Week6_3, file = "./filtered_Week6_3.rds")
saveRDS(Week6_4, file = "./filtered_Week6_4.rds")
saveRDS(Week6_5, file = "./filtered_Week6_5.rds")
saveRDS(Week6_6, file = "./filtered_Week6_6.rds")




#' 
## ----add to metadata of each seurat object------------------------------------------------------------------------------------------------------------
#Assign an identity to the cells - this is only necessary when you are combining more than one file/Seurat object
## Sample Name
#Week 18
Week18_1@meta.data$sample<-"Week18_1"
Week18_2@meta.data$sample<-"Week18_2"
Week18_3@meta.data$sample<-"Week18_3"
Week18_4@meta.data$sample<-"Week18_4"
Week18_5@meta.data$sample<-"Week18_5"
Week18_6@meta.data$sample<-"Week18_6"
#Week 12
Week12_1@meta.data$sample<-"Week12_1"
Week12_2@meta.data$sample<-"Week12_2"
Week12_3@meta.data$sample<-"Week12_3"
Week12_4@meta.data$sample<-"Week12_4"
#Week 6
Week6_1@meta.data$sample<-"Week6_1"
Week6_2@meta.data$sample<-"Week6_2"
Week6_3@meta.data$sample<-"Week6_3"
Week6_4@meta.data$sample<-"Week6_4"
Week6_5@meta.data$sample<-"Week6_5"
Week6_6@meta.data$sample<-"Week6_6"

## Time
#Week 18
Week18_1@meta.data$time<-"Week18"
Week18_2@meta.data$time<-"Week18"
Week18_3@meta.data$time<-"Week18"
Week18_4@meta.data$time<-"Week18"
Week18_5@meta.data$time<-"Week18"
Week18_6@meta.data$time<-"Week18"
#Week 12
Week12_1@meta.data$time<-"Week12"
Week12_2@meta.data$time<-"Week12"
Week12_3@meta.data$time<-"Week12"
Week12_4@meta.data$time<-"Week12"
#Week 6
Week6_1@meta.data$time<-"Week6"
Week6_2@meta.data$time<-"Week6"
Week6_3@meta.data$time<-"Week6"
Week6_4@meta.data$time<-"Week6"
Week6_5@meta.data$time<-"Week6"
Week6_6@meta.data$time<-"Week6"

## Group
#Week 1
Week18_1@meta.data$group<-"Pre-Progressor"
Week18_2@meta.data$group<-"Pre-Progressor"
Week18_3@meta.data$group<-"Pre-Progressor"
Week18_4@meta.data$group<-"Progressor"
Week18_5@meta.data$group<-"Progressor"
Week18_6@meta.data$group<-"Progressor"
#Week 12
Week12_1@meta.data$group<-"Progressor"
Week12_2@meta.data$group<-"Progressor"
Week12_3@meta.data$group<-"Non-Progressor"
Week12_4@meta.data$group<-"Non-Progressor"
#Week 6
Week6_1@meta.data$group<-"Non-Progressor"
Week6_2@meta.data$group<-"Progressor"
Week6_3@meta.data$group<-"Non-Progressor"
Week6_4@meta.data$group<-"Progressor"
Week6_5@meta.data$group<-"Progressor"
Week6_6@meta.data$group<-"Non-Progressor"

#' 
#' # Merge THEN scTransform
#' 
#' Do not run if using scTransform method
#' 
## ----Merge THEN scTransform method--------------------------------------------------------------------------------------------------------------------
#Create Seurat objects
# use count matrix to create seurat object, which has data (ie count matrix)
# and analysis (ie PCA/clustering)
# initialized with raw data

Week18_1 <- NormalizeData(Week18_1) # normalize the data with log norm
Week18_1 <- ScaleData(Week18_1, verbose = F) # linear transformation prior to
# dimensional reduction like PCA. Shifts expression so mean exp is 0 across all cells,
# scales exp, so var across cells is 1
# verbose calls progress bar
Week18_2<- NormalizeData(Week18_2)
Week18_2 <- ScaleData(Week18_2, verbose = F)
Week18_3<- NormalizeData(Week18_3)
Week18_3 <- ScaleData(Week18_3, verbose = F)
Week18_4<- NormalizeData(Week18_4)
Week18_4 <- ScaleData(Week18_4, verbose = F)
Week18_5<- NormalizeData(Week18_5)
Week18_5 <- ScaleData(Week18_5, verbose = F)
Week18_6<- NormalizeData(Week18_6)
Week18_6 <- ScaleData(Week18_6, verbose = F)

Week12_1 <- NormalizeData(Week12_1) 
Week12_1 <- ScaleData(Week12_1, verbose = F) 
Week12_2<- NormalizeData(Week12_2)
Week12_2 <- ScaleData(Week12_2, verbose = F)
Week12_3<- NormalizeData(Week12_3)
Week12_3 <- ScaleData(Week12_3, verbose = F)
Week12_4<- NormalizeData(Week12_4)
Week12_4 <- ScaleData(Week12_4, verbose = F)


Week6_1 <- NormalizeData(Week6_1) # normalize the data with log norm
Week6_1 <- ScaleData(Week6_1, verbose = F) 
Week6_2<- NormalizeData(Week6_2)
Week6_2 <- ScaleData(Week6_2, verbose = F)
Week6_3<- NormalizeData(Week6_3)
Week6_3 <- ScaleData(Week6_3, verbose = F)
Week6_4<- NormalizeData(Week6_4)
Week6_4 <- ScaleData(Week6_4, verbose = F)
Week6_5<- NormalizeData(Week6_5)
Week6_5 <- ScaleData(Week6_5, verbose = F)
Week6_6<- NormalizeData(Week6_6)
Week6_6 <- ScaleData(Week6_6, verbose = F)



#' 
#' merge tutorial: <https://satijalab.org/seurat/archive/v4.3/merge>
#' 
## ----Analyze Data from Week 6 and Week 12--------------------------------------------------------------------------------------------------------------------------------------------
T1D_Timepoints = merge(x = Week6_1, y = list( Week6_2, Week6_3,Week6_4, Week6_5,Week6_6,
                                              Week12_1,Week12_2,Week12_3,Week12_4))


#' 
#' ## scTransform
#' 
#' scTransform vignette:
#' <https://satijalab.org/seurat/articles/sctransform_vignette.html>
#' install glmGamPoi before using, significantly improves speed
#' BiocManager::install("glmGamPoi")
#' 
## ----scTransform on Week 6 and Week 12 timepoints --------------------------------------------------------------------------------------------------------------------------------------
T1D_Timepoints  <- SCTransform(T1D_Timepoints, vars.to.regress = "percent.mt", verbose=F)
saveRDS(T1D_Timepoints , file = "./T1D_Timepoints_sct_filtered.rds")

## ----save merged and cluster unintegrated-------------------------------------------------------------------------------------------------------------
T1D_Timepoints = readRDS("./T1D_Timepoints_sct_filtered.rds")
T1D_Timepoints <- RunPCA(T1D_Timepoints)
T1D_Timepoints <- RunUMAP(T1D_Timepoints, dims = 1:30)
T1D_Timepoints[[]]
DimPlot(T1D_Timepoints, reduction = "umap", group.by = c("sample", "seurat_clusters"))
DimPlot(T1D_Timepoints, reduction = "umap", split.by = c("sample"))

ElbowPlot(T1D_Timepoints, ndims = 50)
T1D_Timepoints <- FindNeighbors(T1D_Timepoints, dims = 1:30, verbose = F)
T1D_Timepoints <- FindClusters(T1D_Timepoints, verbose =F, resolution = 0.4)
DimPlot(T1D_Timepoints, group.by = c("time"), label=F)

DimPlot(T1D_Timepoints, reduction = "umap", group.by = c("sample"))
DimPlot(int.T1D_Timepoints, reduction = "umap", group.by = c("sample"))
table(Idents(T1D_Timepoints),T1D_Timepoints@meta.data$sample)


Layers(T1D_Timepoints)
T1D_Timepoints
DefaultAssay(T1D_Timepoints) = "RNA"
T1D_Timepoints <- JoinLayers(T1D_Timepoints)


get_conserved <- function(cluster){
  Seurat::FindConservedMarkers(T1D_Timepoints,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    cbind(cluster_id = cluster, .)
}
#devtools::install_github('immunogenomics/presto')
FindConservedMarkers(T1D_Timepoints,
                     ident.1 = 0,
                     grouping.var = "sample",
                     only.pos = TRUE)

#conserved_markers <- map_dfr(unique(new.cluster.ids), get_conserved)
DimPlot(T1D_Timepoints,reduction = "umap")
conserved_markers <- map_dfr(0:22, get_conserved)

DefaultAssay(T1D_Timepoints) = "SCT"
head(conserved_markers)
conserved_markers$Week6_2_avg_log2FC
# Extract top 10 markers per cluster
top20 <- conserved_markers %>% 
  mutate(avg_fc = (`Week6_1_avg_log2FC` + 
                     `Week6_2_avg_log2FC` + 
                     `Week6_3_avg_log2FC` +
                     `Week6_4_avg_log2FC` +
                     `Week6_5_avg_log2FC` +
                     `Week6_6_avg_log2FC` +
                     `Week12_1_avg_log2FC` + 
                     `Week12_2_avg_log2FC` +
                     `Week12_3_avg_log2FC` +
                     `Week12_4_avg_log2FC` ) /10) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)




## Remaining Annotation and Clustering Done through mapping with conventional marker and PrepSCTFindMarkers

