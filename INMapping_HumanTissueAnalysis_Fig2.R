# 1. Load Libraries ----
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(biomaRt)
BiocManager::install("glmGamPoi")
library(glmGamPoi)
getwd()


# # Extract raw counts and metadata to create SingleCellExperiment object
# counts_blood <- T1D_Blood@assays$RNA$counts 
# metadata_blood <- T1D_Blood@meta.data
# # Set up metadata as desired for aggregation and DE analysis
# metadata_blood$cluster_id <- factor(T1D_Blood@active.ident)
# # Create single cell experiment object
# sce_blood <- SingleCellExperiment(assays = list(counts = counts_blood), 
#                                colData = metadata_blood)
# dim(colData(sce_blood))
# head(colData(sce_blood))


#2. Map Early-Intermediate DEGs Genes ----

# EarlyStageDEGs<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Early_NOD_ProgressorVsNonProgressor.csv")
# EarlyStageDEGs<-EarlyStageDEGs[abs(EarlyStageDEGs$log2FoldChange)>=1 & EarlyStageDEGs$pvalue<=0.05,]
# 
# IntermediateStageDEGs<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Intermediate_NOD_ProgressorVsNonProgressor.csv")
# IntermediateStageDEGs<-IntermediateStageDEGs[abs(IntermediateStageDEGs$log2FoldChange)>=1 & IntermediateStageDEGs$pvalue<=0.05,]
# 
# LateStageDEGs<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Late_NOD_ProgressorVsNonProgressor.csv")
# LateStageDEGs<-LateStageDEGs[abs(LateStageDEGs$log2FoldChange)>=1 & LateStageDEGs$pvalue<=0.05,]
# 
# early_genes <- unique(EarlyStageDEGs[[1]])
# intermediate_genes <- unique(IntermediateStageDEGs[[1]])
# late_genes <- unique(LateStageDEGs[[1]])
# common_genes_all <- Reduce(intersect, 
#                             list(early_genes, intermediate_genes, late_genes))
# 
# length(common_genes_all)
# head(common_genes_all)
# 
# library(VennDiagram)
# library(grid)
# 
# # Calm, Nature-style palette
# custom_colors <- c(
#   "#4DAC26",  
#   "#C77C57" ,  
#   "#9E4E5E" 
# )
# venn.plot <- venn.diagram(
#   x = list(
#     `Early Stage`        = early_genes,
#     `Intermediate Stage` = intermediate_genes,
#     `Late Stage`         = late_genes
#   ),
#   filename = NULL,
#   fill = custom_colors,
#   alpha = 0.7,
#   lty = "blank",           # No borders
#   cex = 2.4,               # Intersection numbers
#   fontface = "bold",
#   fontfamily = "Helvetica",
#   cat.col = custom_colors,
#   cat.cex = 2.2,
#   cat.fontface = "bold",
#   cat.fontfamily = "Helvetica",
#   margin = 0.08
# )
# 
# grid.newpage()
# grid.draw(venn.plot)
# 
# 
# genes_Early_Intermediate <- intersect(early_genes, intermediate_genes)
# length(genes_Early_Intermediate)
# genes_Early_Intermediate
# 
# 
# genes_Early_Late <- intersect(early_genes, late_genes)
# length(genes_Early_Late)
# genes_Early_Late
# 
# genes_Intermediate_Late <- intersect(intermediate_genes, late_genes)
# length(genes_Intermediate_Late)
# genes_Intermediate_Late
# 
# Persistent_DEGs <- union(
#   genes_Early_Intermediate,
#   genes_Early_Late
# )
# 
# Persistent_DEGs <- union(
#   Persistent_DEGs,
#   genes_Intermediate_Late
# )
# 
# 
# length(Persistent_DEGs)
# Persistent_DEGs
# 
# Persistent_Immune_DEGs <- c(
#   "Cxcl15",   # chemokine (neutrophil recruitment)
#   "Apol6",    # interferon-inducible, apoptosis/innate immunity
#   "Itih4",    # acute phase protein
#   "Wfdc18",   # immune / mucosal defense-associated
#   "Ccr8",     # T cell / Treg chemokine receptor
#   "Efna3",    # immune cell migration & immune–stromal signaling
#   "Fgg",      # coagulation–immunity cross-talk
#   "Vtcn1",    # B7-H4 immune checkpoint
#   "Agtr2",    # immunomodulatory (macrophage, T cell signaling)
#   "C6",       # complement
#   "Cd7",      # T cell / NK marker
#   "Klrd1",    # NK cell receptor (CD94)
#   "Ptgdr",    # prostaglandin D2 receptor (Th2, DCs)
#   "Rorc"      # Th17 / ILC3 master transcription factor
# )


#Select Immune and Inflammtion DEGs from Early and Intermediate Timepoiny
DEG_immune_genes_mouse <- c(
  ### Early Stage
  # Neutrophils / innate
   "S100a8", "S100a9",
  # Myeloid / macrophage
  "Cd163l1", "Itih4",
  # B cells (adaptive)
  "Ms4a1", "Blk",
  # T cells (adaptive)
  "Cd6", "Cd40lg", "Ctla4", "Ccr8", "Il23r", "Themis",
  ### Intermediate Stage
  # Myeloid / inflammatory
  "Cd5l", "Lilra5", "Fpr1", "C6",
  # Interferon / ISGs
  "Ifit1", "Ifit2", "Ifit3", "Oas3", "Rtp4",
  # NK / cytotoxic
  "Klrd1", "Xcl1",
  # B cells / antibody
  "Ighg1",
  # T cells / lymphoid (adaptive – first)
  "Cd4", "Cd69", "Cd7", "Gzmk", "Rorc"
)

DEG_immune_genes_human <- c(
  ### Early Stage
  # Neutrophil / innate
  "S100A8", "S100A9",
  # Myeloid / macrophage
  "CD163L1", "ITIH4",
  # B cells (adaptive)
  "MS4A1", "BLK",
  # T cells (adaptive)
  "CD6", "CD40LG", "CTLA4", "CCR8", "IL23R", "THEMIS",
  ### Intermediate Stage
  # Myeloid / innate
  "CD5L", "LILRA5", "FPR1", "C6",
  # Interferon / ISGs
  "IFIT1", "IFIT2", "IFIT3", "OAS3", "RTP4",
  # NK / cytotoxic
  "KLRD1", "XCL1",
  # B cells
  "IGHG1",
  # T cells / lymphoid
  "CD4", "CD69", "CD7", "GZMK", "RORC"
)


## NOD Pancreas Immune  Cells----
NOD_T1D_Timepoints<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v8.rds")

tcell_idents <- c(
  "CD8 memory",
  "CD8 exhausted effector-like",
  "Tcon memory",
  "Tcon activated ",
  "Tcon exhausted effector-like",
  "Tcon Interferon Sensing",
  "Th17-like",
  "Th2-like",
  "Tregs",
  "Gamma Delta T Cell",
  "Cytotoxic NK/CD8"
)

new_ids <- as.character(Idents(NOD_T1D_Timepoints))
new_ids[new_ids %in% tcell_idents] <- "T Cell"
Idents(NOD_T1D_Timepoints) <- factor(new_ids)
unique(Idents(NOD_T1D_Timepoints))

Idents(NOD_T1D_Timepoints) <- factor(
  Idents(NOD_T1D_Timepoints),
  levels = c(
    "B Cell",
    "T Cell",
    "Acinar Cell",
    "Dendritic Cell",
    "Macrophage",
    "Plasma Cell"
  )
)

levels(Idents(NOD_T1D_Timepoints))


DotPlot(NOD_T1D_Timepoints, features = DEG_immune_genes_mouse, dot.scale = 10) +
  labs(
    y = "Pancreas Cell Type",
    x = "Niche DEGs"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )+
  scale_colour_gradient2(
    low = "navyblue",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  )+
  scale_size(
    limits = c(0, 100),
    range = c(0, 10)
  ) 



## Blood Immune Cells----

T1D_Blood = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/Blood/T1D_Seurat_Object_Final_v2.rds")
T1D_Blood <- subset(T1D_Blood, idents = setdiff(levels(Idents(T1D_Blood)), "Undetermined"))

unique(Idents(T1D_Blood))

T1D_Blood <- RenameIdents(
  T1D_Blood,
  "B-Plasma"="Plasma Cell",
  "NK"="NK Cell",
  "cDCs"="Dendritic Cell",
  "pDCs"="Dendritic Cell",
  "B-Naive"="B Cell",
  "B-SM"="B Cell",
  "T-reg"="T Cell",
  "CD4-CM"="T Cell",
  "CD4-EM"="T Cell",
  "CD8-EM"="T Cell",
  "CD8-CM"="T Cell",
  "CD4-Naive"="T Cell",
  "CD8-Naive"="T Cell",
  "VD2p"="T Cell",
  "MAIT"="T Cell",
  "C-Monocyte"="Monocyte",
  "I-Monocyte"="Monocyte",
  "NC-Monocyte"="Monocyte"
)

unique(Idents(T1D_Blood))
Idents(T1D_Blood) <- factor(
  Idents(T1D_Blood),
  levels = c(
    "B Cell",
    "T Cell",
    "HSCs",
    "Erythrocytes",
    "Platelet",
    "NK Cell",
    "Dendritic Cell",
    "Monocyte",
    "Plasma Cell"
  )
)



DimPlot(T1D_Blood,label = TRUE, label.box = T,label.size = 9,repel = T)+
  NoAxes() +NoLegend()


DotPlot(T1D_Blood, features = DEG_immune_genes_human, dot.scale = 10) +
  labs(
    y = "Blood Cell Type",
    x = "Niche DEGs (Human Orthologs)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )+
  scale_colour_gradient2(
    low = "navyblue",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_size(
    limits = c(0, 100),
    range = c(0, 10)
  ) 

## Spleen  Cells----
T1D_Spleen = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/MAI_T1D/spleen_result_annotated_v2.rds")
T1D_Spleen <- RenameIdents(
  T1D_Spleen,
  "B" = "B Cell",
  "T" = "T Cell",
  "Cycling" = "Cycling Cell",
  "NK"="NK Cell"
)

DimPlot(T1D_Spleen,label = TRUE, label.box = T,label.size = 9,repel = T)+
  NoAxes() +NoLegend()
DEG_immune_genes_human_spleen<- DEG_immune_genes_human[
  !is.na(DEG_immune_genes_human) &
    DEG_immune_genes_human %in% rownames(T1D_Spleen)
]
DotPlot(T1D_Spleen, features = DEG_immune_genes_human_spleen, dot.scale = 10) +
  labs(
    y = "Spleen Cell Type",
    x = "Niche DEGs (Human Orthologs)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )+
  scale_colour_gradient2(
    low = "navyblue",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_size(
    limits = c(0, 100),
    range = c(0, 10)
  ) 



## plN  Cells----

T1D_pln = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/MAI_T1D/pln_result_annotated_v2.rds")
T1D_pln <- RenameIdents(
  T1D_pln,
  "Bcell" = "B Cell",
  "Tcell" = "T Cell",
  "Cycling" = "Cycling Cell",
  "NK"="NK Cell"
)

DimPlot(T1D_pln,label = TRUE, label.box = T,label.size = 9,repel = T)+
  NoAxes()+NoLegend()


DotPlot(T1D_pln, features = DEG_immune_genes_human, dot.scale = 10) +
  labs(
    y = "pLN Cell Type",
    x = "Niche DEGs (Human Orthologs)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )+
  scale_colour_gradient2(
    low = "navyblue",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_size(
    limits = c(0, 100),
    range = c(0, 10)
  ) 





## Islet  Cells----
options(future.globals.maxSize = 20 * 1024^3)library(future)
library(future)
plan(sequential)  # avoids extra memory duplication

T1D_Islet = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/MAI_T1D/islet_final_annotated.rds")
DimPlot(T1D_Islet,label = TRUE, label.box = T,label.size = 6,repel = T)+
  NoAxes() 
T1D_Islet <- SCTransform(T1D_Islet,conserve.memory = T, verbose = T)
saveRDS(T1D_Islet,"/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/MAI_T1D/islet_final_annotated_v2.rds")
DefaultAssay(T1D_Islet)<-"RNA"
DotPlot(T1D_Islet, features = human_genes, dot.scale = 10) +
  scale_colour_gradient2(
    low = "navyblue",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0
  ) +
  coord_flip() +
  scale_x_discrete(labels = function(x) ifelse(x %in% genes_to_label_human, x, "")) +
  labs(
    y = "Blood Cell Type",
    x = "Gene Signature (Human Orthologs)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
    axis.text.y = element_text(size = 10)
  )


#2. Map Top 100 Genes ----
T1DGene_top100_humanorthologs<- read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Input/SVC Model/Top100Genes_T1D_HumanOrthologs.csv", sep=",", header=T,check.names = FALSE)
Top100_Mouse<- T1DGene_top100_humanorthologs$mouse_symbol
Top100_Human<- T1DGene_top100_humanorthologs$human_symbol

## NOD Pancreas Immune  Cells----
NOD_T1D_Timepoints<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v8.rds")

genes_to_label <- c(
  "Pdcd1", "Cd40lg", "Ccr8", "Cpvl", "Cfi", "Tnfaip8l1",
  "Cpt2", "Gpd2", "Aacs", "Mpst",
  "Ncor1", "Sp8", "Yy2", "Zfp703",
  "Nedd4l", "Psme3", "Rab3c", "Snx3",
  "Col10a1", "Bambi",
  "Zcchc18",
  "Tm9sf4", "Cdkn2c", "Stmn1", "Lbh", "Rbp4",
  "Dusp12", "Prkcz2", "Gsk3a", "Wnk1",
  "Cth",
  "Ephx3", "Pgrmc2", "Nppa",
  "Ywhag", #Mac  signaling
  "Apol10b",#Endothelial marker present in macrophage
  "B3gnt5",#B cells
  "Dynll2"# Involved In DC intracellular transport and CD8 T cell ant-tumor immunity
  )

Top100_Mouse_filtered <- Top100_Mouse[
  !is.na(Top100_Mouse) &
    Top100_Mouse %in% rownames(NOD_T1D_Timepoints)
]
#65 genes

DotPlot(NOD_T1D_Timepoints, features = Top100_Mouse_filtered, dot.scale = 10) +
  labs(
    y = "NOD Pancreas",
    x = "IN Gene Signature"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )+
  scale_colour_gradient2(
    low = "navyblue",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_size(
    limits = c(0, 100),
    range = c(0, 10)
  ) +coord_flip()+
  scale_x_discrete(
    labels = function(y) ifelse(y %in% genes_to_label, y, "")
  ) 

## Blood Immune Cells----

#T1D_Blood = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/Blood/T1D_Seurat_Object_Final_v2.rds")

Top100_Mouse<- T1DGene_top100_humanorthologs$mouse_symbol
genes_label_human <- T1DGene_top100_humanorthologs[T1DGene_top100_humanorthologs$mouse_symbol %in% genes_to_label,]$human_symbol
Top100_Human_blood<- Top100_Human[
  !is.na(Top100_Human) &
    Top100_Human %in% rownames(T1D_Blood)
]
#73 genes
DotPlot(T1D_Blood, features = Top100_Human_blood, dot.scale = 10) +
  labs(
    y = "Human Blood",
    x = "IN Gene Signature"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )+
  scale_colour_gradient2(
    low = "navyblue",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_size(
    limits = c(0, 100),
    range = c(0, 10)
  ) +coord_flip()+
  scale_x_discrete(
    labels = function(y) ifelse(y %in% genes_label_human, y, "")
  ) 

## Spleen  Cells----

#T1D_Spleen = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/MAI_T1D/spleen_result_annotated_v2.rds")
Top100_Human_spleen<- Top100_Human[
  !is.na(Top100_Human) &
    Top100_Human %in% rownames(T1D_Spleen)
]
#89 genes
DotPlot(T1D_Spleen, features = Top100_Human_blood, dot.scale = 10) +
  labs(
    y = "Human Spleen",
    x = "IN Gene Signature"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )+
  scale_colour_gradient2(
    low = "navyblue",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_size(
    limits = c(0, 100),
    range = c(0, 10)
  ) +coord_flip()+
  scale_x_discrete(
    labels = function(y) ifelse(y %in% genes_label_human, y, "")
  ) 

## plN  Cells----


T1D_pln = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/MAI_T1D/pln_result_annotated_v2.rds")

Top100_Human_pln<- Top100_Human[
  !is.na(Top100_Human) &
    Top100_Human %in% rownames(T1D_pln)
]
#87 genes
DotPlot(T1D_pln, features = Top100_Human_pln, dot.scale = 10) +
  labs(
    y = "Human pLN",
    x = "IN Gene Signature"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  )+
  scale_colour_gradient2(
    low = "navyblue",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  scale_size(
    limits = c(0, 100),
    range = c(0, 10)
  ) +coord_flip()+
  scale_x_discrete(
    labels = function(y) ifelse(y %in% genes_label_human, y, "")
  ) 


