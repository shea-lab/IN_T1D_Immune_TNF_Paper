## ----load libraries-----------------------------------------------------------------------------------------------------------------------------------
setwd("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/SingleCellRNASeq/ProcessedRDSFile")
library(Seurat)
packageVersion("Seurat")
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
library(glmGamPoi)
library(msigdbr)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(forcats)
library(EnhancedVolcano)
library(DESeq2)
library(Matrix.utils)
library(enrichplot)

suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(SingleCellExperiment))
#BiocManager::install("scran")
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratObject))
library(clusterProfiler)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)  # Use org.Hs.eg.db for human genes
library(org.Hs.eg.db) 

# Hallmark
hallmark <- msigdbr(species = "Mus musculus", category  = "H")
mm_h_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
mm_h_df <- data.frame(
  gs_name = rep(names(mm_h_sets), sapply(mm_h_sets, length)),
  gene_symbol = unlist(mm_h_sets)
)


# KEGG
kegg_all <- msigdbr(species="Mus musculus", category="C2", subcategory="CP:KEGG_LEGACY")
mm_kegg_sets <- split(kegg_all$gene_symbol, kegg_all$gs_name)
mm_kegg_df <- data.frame(
  gs_name = rep(names(mm_kegg_sets), sapply(mm_kegg_sets, length)),
  gene_symbol = unlist(mm_kegg_sets)
)

mm_all_df <- rbind(mm_h_df, mm_kegg_df)

## Human version
hallmark_h <- msigdbr(species = "Homo sapiens", category  = "H")
mm_h_sets_h <- split(hallmark_h$gene_symbol, hallmark_h$gs_name)
mm_h_df_h <- data.frame(
  gs_name = rep(names(mm_h_sets_h), sapply(mm_h_sets_h, length)),
  gene_symbol = unlist(mm_h_sets_h)
)

kegg_all_h <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:KEGG_LEGACY")
mm_kegg_sets_h <- split(kegg_all_h$gene_symbol, kegg_all_h$gs_name)
mm_kegg_df_h <- data.frame(
  gs_name = rep(names(mm_kegg_sets_h), sapply(mm_kegg_sets_h, length)),
  gene_symbol = unlist(mm_kegg_sets_h)
)
mm_all_df_human <- rbind(mm_h_df_h, mm_kegg_df_h)



## NOD Pancreas----
NOD_T1D_Timepoints<-readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/ProcessedRDSFile/annotated_T1D_Timepoints_v8.rds")
NOD_T1D_Timepoints
Macrophage_NOD <- subset(NOD_T1D_Timepoints, subset = (CellSubType == "Macrophage"))
DefaultAssay(Macrophage_NOD) <- "RNA"
counts_mac_NOD <- GetAssayData(Macrophage_NOD, assay = "RNA", layer = "counts")
sample_id_mac_NOD <- factor(Macrophage_NOD$sample)  

pb_mac_NOD <- aggregate.Matrix(
  t(counts_mac_NOD),
  groupings = data.frame(sample = sample_id_mac_NOD),
  fun = "sum"
)
pb_mac_NOD<-t(pb_mac_NOD)

#### DESEQ----

meta_mac_NOD <- data.frame(
  row.names = colnames(pb_mac_NOD),
  group = Macrophage_NOD@meta.data$group[
    match(colnames(pb_mac_NOD), Macrophage_NOD$sample)
  ],
  stage = Macrophage_NOD@meta.data$time[
    match(colnames(pb_mac_NOD), Macrophage_NOD$sample)
  ]
)

meta_mac_NOD$group <- factor(meta_mac_NOD$group,
                             levels = c("Non-Progressor", "Progressor"))

meta_mac_NOD$stage <- factor(meta_mac_NOD$stage,
                             levels = c("Week6", "Week12"))



#Remove all zero genes
pb_mac_NOD <- pb_mac_NOD[rowSums(pb_mac_NOD) > 0, ]

dds_mac_NOD <- DESeqDataSetFromMatrix(
  countData = round(pb_mac_NOD),
  colData   = meta_mac_NOD,
  design    = ~ stage+ group
)

dds_mac_NOD <- DESeq(dds_mac_NOD)
res_mac_NOD <- results(dds_mac_NOD)

res_mac_NOD_df <- as.data.frame(res_mac_NOD)
res_mac_NOD_df$pvalue[is.na(res_mac_NOD_df$pvalue)] <- 1
res_mac_NOD_df$padj[is.na(res_mac_NOD_df$padj)] <- 1
deg_genes_mac_NOD <- rownames(res_mac_NOD_df)[
  res_mac_NOD_df$pvalue < 0.05 &
    abs(res_mac_NOD_df$log2FoldChange) >= 1
]

# Define color mapping with stronger contrast for visibility
keyvals <- ifelse(res_mac_NOD_df$log2FoldChange > 1 & res_mac_NOD_df$pvalue < 0.05, "#FF0000",  # Bright red for nominally significant up
                  ifelse(res_mac_NOD_df$log2FoldChange < -1 & res_mac_NOD_df$pvalue < 0.05, "#1E90FF",  # Bright blue for nominally significant down
                         "gray"))  # Gray for non-significant genes

# Assign names for legend
names(keyvals) <- ifelse(res_mac_NOD_df$log2FoldChange > 1 & res_mac_NOD_df$pvalue < 0.05, "p-value < 0.05 & log[2]FC > 1",
                         ifelse(res_mac_NOD_df$log2FoldChange < -1 & res_mac_NOD_df$pvalue < 0.05, "p-value < 0.05 & log[2]FC < -1",
                                "Not Significant"))

point_sizes <-  2.0


# Genes to label
genes_to_label <- c(
  "Tnfaip3",
  "Ccl5",
  "Dusp1",
  "Pde4b",
  "C5ar2",
  "Cd83",
  "Mertk",
  "Il6ra",
  "Irs2",
  "Il1rn",
  "Anxa1",
  "Cd200",
  "Socs6",
  "Cryab",
  "Selenos",
  "Clu",
  "Aldh1a2"
)


EnhancedVolcano(res_mac_NOD_df,
                lab = rownames(res_mac_NOD_df),
                selectLab = genes_to_label,
                x = 'log2FoldChange',
                y = 'pvalue',  
                title = 'NOD Pancreas Macrophage',
                pCutoff = 0.05,  
                FCcutoff = 1, 
                pointSize = point_sizes,  
                labSize = 8,
                colAlpha = 0.75,
                legendLabels = c("Not Significant", "p-value < 0.05"),
                legendPosition = 'right',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                boxedLabels = TRUE,
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  
                ylab = expression("-log"[10] ~ "(p-value)"),  
                xlim = c(-8, 8), 
                ylim = c(0, 6),  
                max.overlaps = 40,
                colCustom = keyvals, 
                caption = NULL 
)



# Save the data frame as CSV
write.csv(res_mac_NOD_df, file = "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 4/Data-Results/DEG_DESEQ_NODPancreasMacrophage_ProgressorVsNonProgressor.csv")

#### GSEA----

#Build ranked gene lists using DESeq2 Wald stat
lfc_vector_mac_NOD  <- res_mac_NOD_df$stat;  names(lfc_vector_mac_NOD)  <- rownames(res_mac_NOD_df)

# Drop NAs
lfc_vector_mac_NOD  <- lfc_vector_mac_NOD[!is.na(lfc_vector_mac_NOD)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_mac_NOD  <- sort(lfc_vector_mac_NOD,  decreasing = TRUE)


gsea_results_mac_NOD <- GSEA(
  geneList = lfc_vector_mac_NOD, 
  minGSSize = 5, 
  maxGSSize = 500, # Maximum gene set size
  pvalueCutoff = 1, 
  eps = 0, 
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = mm_all_df  
)

# Extract results for NOD Pacreas Macrophages
gsea_results_mac_NOD_df <- as.data.frame(gsea_results_mac_NOD)

write.csv(gsea_results_mac_NOD_df, file = "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 4/Data-Results/GSEA_NODPancreasMacrophage_ProgressorVsNonProgressor.csv")

library(enrichplot)

gseaplot2(
  gsea_results_mac_NOD,
  geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  title = "TNFα Signaling via NF-KB",
  base_size = 25
)

tnf_row <- gsea_results_mac_NOD_df[
  gsea_results_mac_NOD_df$Description == "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
]

# Human Islets----
T1D_Islets = readRDS( "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/MAI_T1D/hpap_islet.rds")
unique(T1D_Islets@meta.data$`Diabetes Status`)
T1D_Islets <- subset(
  T1D_Islets,
  subset = `Diabetes Status` != "T2D"
)
unique(T1D_Islets@meta.data$`Diabetes Status`)
Idents(T1D_Islets)<-T1D_Islets$`Cell Type Grouped`
DimPlot(T1D_Islets,label = T, pt.size = 1.2,label.box = T,repel = T,label.size = 7)+
  NoLegend() +   
  NoAxes() +       
  theme(
    axis.line = element_blank()
  )

DimPlot(T1D_Islets, pt.size = 1.2,group.by ="Diabetes Status" )+      # removes legend
  NoAxes() +       
  theme(
    axis.line = element_blank()
  )


#### Pseudobulk----
Macrophage_Islet <- subset(T1D_Islets, subset = (`Cell Type Grouped` == "Macrophage"))
DefaultAssay(Macrophage_Islet) <- "RNA"
counts_Mac_Islets <- GetAssayData(Macrophage_Islet, assay = "RNA", layer = "counts")
sample_id_Mac_Islets <- factor(Macrophage_Islet$Library)  

pb_mac_Islet <- aggregate.Matrix(
  t(counts_Mac_Islets),
  groupings = data.frame(sample = sample_id_Mac_Islets),
  fun = "sum"
)
pb_mac_Islet<-t(pb_mac_Islet)
head(pb_mac_Islet)

#### DESEQ----

meta_mac_Islet <- data.frame(
  row.names = colnames(pb_mac_Islet),
  group = Macrophage_Islet@meta.data$`Diabetes Status`[
    match(colnames(pb_mac_Islet), Macrophage_Islet$Library)
  ]
)

meta_mac_Islet$group <- factor(meta_mac_Islet$group,
                                levels = c("ND","AAB+","T1D"))

#Remove all zero genes
pb_mac_Islet <- pb_mac_Islet[rowSums(pb_mac_Islet) > 0, ]

dds_mac_Islet <- DESeqDataSetFromMatrix(
  countData = round(pb_mac_Islet),
  colData   = meta_mac_Islet,
  design    = ~  group
)

dds_mac_Islet <- DESeq(dds_mac_Islet)
results(dds_mac_Islet)
# 1) T1D vs ND
res_T1D_vs_ND <- results(dds_mac_Islet, contrast = c("group", "T1D", "ND"))

# 2) T1D vs AAB+
res_T1D_vs_AAB <- results(dds_mac_Islet, contrast = c("group", "T1D", "AAB+"))

res_T1D_vs_ND_df <- as.data.frame(res_T1D_vs_ND)
res_T1D_vs_ND_df$pvalue[is.na(res_T1D_vs_ND_df$pvalue)] <- 1
res_T1D_vs_ND_df$padj[is.na(res_T1D_vs_ND_df$padj)] <- 1


res_T1D_vs_AAB_df <- as.data.frame(res_T1D_vs_AAB)
res_T1D_vs_AAB_df$pvalue[is.na(res_T1D_vs_AAB_df$pvalue)] <- 1
res_T1D_vs_AAB_df$padj[is.na(res_T1D_vs_AAB_df$padj)] <- 1


# Save the data frame as CSV
write.csv(res_T1D_vs_ND_df, file = "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 4/Data-Results/DEG_DESEQ_HumanMacrophageIslets_T1DVsND.csv")
write.csv(res_T1D_vs_AAB_df, file = "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 4/Data-Results/DEG_DESEQ_HumanMacrophageIslets_T1DVsAAB.csv")

#### GSEA----

##### T1D vs ND----
#Build ranked gene lists using DESeq2 Wald stat
lfc_vector_T1D_vs_ND  <- res_T1D_vs_ND_df$stat;  names(lfc_vector_T1D_vs_ND)  <- rownames(res_T1D_vs_ND_df)

# Drop NAs
lfc_vector_T1D_vs_ND  <- lfc_vector_T1D_vs_ND[!is.na(lfc_vector_T1D_vs_ND)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_T1D_vs_ND  <- sort(lfc_vector_T1D_vs_ND,  decreasing = TRUE)


gsea_results_mac_islet_T1D_vs_ND <- GSEA(
  geneList = lfc_vector_T1D_vs_ND, 
  minGSSize = 5, 
  maxGSSize = 500, 
  pvalueCutoff = 1, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = mm_all_df_human 
)

# Extract results for Huamn Islet Macrophages - T1D vs  ND
gsea_results_mac_islet_T1D_vs_ND_df <- as.data.frame(gsea_results_mac_islet_T1D_vs_ND)

write.csv(gsea_results_mac_islet_T1D_vs_ND_df, file = "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 4/Data-Results/GSEA_HumanMacrophageIslets_T1DVsND.csv")


gseaplot2(
  gsea_results_mac_islet_T1D_vs_ND,
  geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  title = "TNFα Signaling via NF-KB",
  base_size = 25
)

tnf_row <- gsea_results_mac_islet_T1D_vs_ND[
  gsea_results_mac_islet_T1D_vs_ND$Description == "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
]

##### T1D vs AAB----
#Build ranked gene lists using DESeq2 Wald stat
lfc_vector_T1D_vs_AAB  <- res_T1D_vs_AAB_df$stat;  names(lfc_vector_T1D_vs_AAB)  <- rownames(res_T1D_vs_AAB_df)

# Drop NAs
lfc_vector_T1D_vs_AAB  <- lfc_vector_T1D_vs_AAB[!is.na(lfc_vector_T1D_vs_AAB)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_T1D_vs_AAB  <- sort(lfc_vector_T1D_vs_AAB,  decreasing = TRUE)


gsea_results_mac_islet_T1D_vs_AAB <- GSEA(
  geneList = lfc_vector_T1D_vs_AAB, 
  minGSSize = 5,
  maxGSSize = 500, 
  pvalueCutoff = 1,
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = mm_all_df_human 
)

# Extract results for Huamn Islet Macrophages - T1D vs AAB
gsea_results_mac_islet_T1D_vs_AAB_df <- as.data.frame(gsea_results_mac_islet_T1D_vs_AAB)

write.csv(gsea_results_mac_islet_T1D_vs_AAB_df, file = "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 4/Data-Results/GSEA_HumanMacrophageIslets_T1DVsAAB.csv")

gseaplot2(
  gsea_results_mac_islet_T1D_vs_AAB,
  geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  title = "TNFα Signaling via NF-KB",
  base_size = 25
)

tnf_row_AAB <- gsea_results_mac_islet_T1D_vs_AAB[
  gsea_results_mac_islet_T1D_vs_AAB$Description == "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
]





