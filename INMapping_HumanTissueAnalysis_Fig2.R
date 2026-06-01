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


# 2. Map DEGs ----

#Select Immune and Inflammtion IN based DEGs from Early and Intermediate Stages
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
DimPlot(NOD_T1D_Timepoints)
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


# 3. DEGs Module Scores----

## NOD Pancreas----

EarlyDEG_NOD<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Early_NOD_ProgressorVsNonProgressor.csv",row.names = 1)
Early_upDEG_NOD <-  EarlyDEG_NOD[EarlyDEG_NOD$log2FoldChange>=1 & EarlyDEG_NOD$pvalue<=0.05, ]
Early_downDEG_NOD <-  EarlyDEG_NOD[EarlyDEG_NOD$log2FoldChange<=-1 & EarlyDEG_NOD$pvalue<=0.05, ]
Early_upDEG_NOD_use <- intersect(rownames(Early_upDEG_NOD), rownames(NOD_T1D_Timepoints))
Early_downDEG_NOD_use <- intersect(rownames(Early_downDEG_NOD), rownames(NOD_T1D_Timepoints))


IntermediateDEG_NOD<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Intermediate_NOD_ProgressorVsNonProgressor.csv",row.names = 1)
Intermediate_upDEG_NOD <-  IntermediateDEG_NOD[IntermediateDEG_NOD$log2FoldChange>=1 & IntermediateDEG_NOD$pvalue<=0.05, ]
Intermediate_downDEG_NOD <-  IntermediateDEG_NOD[IntermediateDEG_NOD$log2FoldChange<=-1 & IntermediateDEG_NOD$pvalue<=0.05, ]
Intermediate_upDEG_NOD_use <- intersect(rownames(Intermediate_upDEG_NOD), rownames(NOD_T1D_Timepoints))
Intermediate_downDEG_NOD_use <- intersect(rownames(Intermediate_downDEG_NOD), rownames(NOD_T1D_Timepoints))

LateDEG_NOD<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Late_NOD_ProgressorVsNonProgressor.csv",row.names = 1)
Late_upDEG_NOD <-  LateDEG_NOD[LateDEG_NOD$log2FoldChange>=1 & LateDEG_NOD$pvalue<=0.05, ]
Late_downDEG_NOD <-  LateDEG_NOD[LateDEG_NOD$log2FoldChange<=-1 & LateDEG_NOD$pvalue<=0.05, ]
Late_upDEG_NOD_use <- intersect(rownames(Late_upDEG_NOD), rownames(NOD_T1D_Timepoints))
Late_downDEG_NOD_use <- intersect(rownames(Late_downDEG_NOD), rownames(NOD_T1D_Timepoints))


# Add module scores

#Early
NOD_T1D_Timepoints <- AddModuleScore(
  object = NOD_T1D_Timepoints,
  features = list(Early_upDEG_NOD_use),
  name = "EarlyProgressor_Score",
  assay = DefaultAssay(NOD_T1D_Timepoints)
)
NOD_T1D_Timepoints <- AddModuleScore(
  object = NOD_T1D_Timepoints,
  features = list(Early_downDEG_NOD_use),
  name = "EarlyNonProgressor_Score",
  assay = DefaultAssay(NOD_T1D_Timepoints)
)

#Intermediate
NOD_T1D_Timepoints <- AddModuleScore(
  object = NOD_T1D_Timepoints,
  features = list(Intermediate_upDEG_NOD_use),
  name = "IntermediateProgressor_Score",
  assay = DefaultAssay(NOD_T1D_Timepoints)
)
NOD_T1D_Timepoints <- AddModuleScore(
  object = NOD_T1D_Timepoints,
  features = list(Intermediate_downDEG_NOD_use),
  name = "IntermediateNonProgressor_Score",
  assay = DefaultAssay(NOD_T1D_Timepoints)
)

#Late
NOD_T1D_Timepoints <- AddModuleScore(
  object = NOD_T1D_Timepoints,
  features = list(Late_upDEG_NOD_use),
  name = "LateProgressor_Score",
  assay = DefaultAssay(NOD_T1D_Timepoints)
)
NOD_T1D_Timepoints <- AddModuleScore(
  object = NOD_T1D_Timepoints,
  features = list(Late_downDEG_NOD_use),
  name = "LateNonProgressor_Score",
  assay = DefaultAssay(NOD_T1D_Timepoints)
)
NOD_T1D_Timepoints$celltype<-Idents(NOD_T1D_Timepoints)


#Early
plot_df <- NOD_T1D_Timepoints@meta.data %>%
  dplyr::select(celltype, EarlyProgressor_Score1, EarlyNonProgressor_Score1) %>%
  filter(!is.na(celltype)) %>%
  pivot_longer(
    cols = c(EarlyProgressor_Score1, EarlyNonProgressor_Score1),
    names_to = "signature",
    values_to = "score"
  )
stat_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(score ~ signature)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = max(score, na.rm = TRUE) + 0.08 * diff(range(plot_df$score, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = paste0(signif(p_adj, 3))
  )


ggplot(plot_df, aes(x = celltype, y = score, fill = signature)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_text(
    data = stat_df,
    aes(x = celltype, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 4.5
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Module Score",
    title = "Early Stage DEGs-NOD Pancreas"
  ) +
  theme(
    axis.title.x = element_text(size = 18),   # x-axis label
    axis.title.y = element_text(size = 18),   # y-axis label
    axis.text.x  = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 18),
    legend.title = element_blank()
  )+scale_fill_manual(
    values = c(
      "EarlyProgressor_Score1" = "#E60000",
      "EarlyNonProgressor_Score1" = "#3D7A60"
    ),
    labels = c(
      "EarlyProgressor_Score1" = "Progressor Score",
      "EarlyNonProgressor_Score1" = "Non-Progressor Score"
    )
  )

median_diff_df <- plot_df %>%
  group_by(celltype, signature) %>%
  summarise(
    median_score = median(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = signature,
    values_from = median_score
  ) %>%
  mutate(
    median_difference = EarlyProgressor_Score1 - EarlyNonProgressor_Score1
  ) %>%
  arrange(desc(abs(median_difference)))


effect_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    median_prog = median(score[signature=="EarlyProgressor_Score1"]),
    median_np = median(score[signature=="EarlyNonProgressor_Score1"]),
    median_diff = median_prog - median_np,
    mad_all = mad(score),
    standardized_effect = median_diff/mad_all
  )


#Intermediate
plot_df <- NOD_T1D_Timepoints@meta.data %>%
  dplyr::select(celltype, IntermediateProgressor_Score1, IntermediateNonProgressor_Score1) %>%
  filter(!is.na(celltype)) %>%
  pivot_longer(
    cols = c(IntermediateProgressor_Score1, IntermediateNonProgressor_Score1),
    names_to = "signature",
    values_to = "score"
  )
stat_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(score ~ signature)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = max(score, na.rm = TRUE) + 0.08 * diff(range(plot_df$score, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = paste0(signif(p_adj, 3))
  )


ggplot(plot_df, aes(x = celltype, y = score, fill = signature)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_text(
    data = stat_df,
    aes(x = celltype, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 4.5
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Module Score",
    title = "Intermediate Stage DEGs-NOD Pancreas"
  ) +
  theme(
    axis.title.x = element_text(size = 18),   # x-axis label
    axis.title.y = element_text(size = 18),   # y-axis label
    axis.text.x  = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 18),
    legend.title = element_blank()
  )+scale_fill_manual(
    values = c(
      "IntermediateProgressor_Score1" = "#E60000",
      "IntermediateNonProgressor_Score1" = "#3D7A60"
    ),
    labels = c(
      "IntermediateProgressor_Score1" = "Progressor Score",
      "IntermediateNonProgressor_Score1" = "Non-Progressor Score"
    )
  )

median_diff_df <- plot_df %>%
  group_by(celltype, signature) %>%
  summarise(
    median_score = median(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = signature,
    values_from = median_score
  ) %>%
  mutate(
    median_difference = IntermediateProgressor_Score1 - IntermediateNonProgressor_Score1
  ) %>%
  arrange(desc(abs(median_difference)))


effect_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    median_prog = median(score[signature=="IntermediateProgressor_Score1"]),
    median_np = median(score[signature=="IntermediateNonProgressor_Score1"]),
    median_diff = median_prog - median_np,
    mad_all = mad(score),
    standardized_effect = median_diff/mad_all
  )


#Late
plot_df <- NOD_T1D_Timepoints@meta.data %>%
  dplyr::select(celltype, LateProgressor_Score1, LateNonProgressor_Score1) %>%
  filter(!is.na(celltype)) %>%
  pivot_longer(
    cols = c(LateProgressor_Score1, LateNonProgressor_Score1),
    names_to = "signature",
    values_to = "score"
  )
stat_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(score ~ signature)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = max(score, na.rm = TRUE) + 0.08 * diff(range(plot_df$score, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    label = paste0(signif(p_adj, 3))
  )


ggplot(plot_df, aes(x = celltype, y = score, fill = signature)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_text(
    data = stat_df,
    aes(x = celltype, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 4.5
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Module Score",
    title = "Late Stage DEGs-NOD Pancreas"
  ) +
  theme(
    axis.title.x = element_text(size = 18),   # x-axis label
    axis.title.y = element_text(size = 18),   # y-axis label
    axis.text.x  = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 18),
    legend.title = element_blank()
  )+scale_fill_manual(
    values = c(
      "LateProgressor_Score1" = "#E60000",
      "LateNonProgressor_Score1" = "#3D7A60"
    ),
    labels = c(
      "LateProgressor_Score1" = "Progressor Score",
      "LateNonProgressor_Score1" = "Non-Progressor Score"
    )
  )

median_diff_df <- plot_df %>%
  group_by(celltype, signature) %>%
  summarise(
    median_score = median(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = signature,
    values_from = median_score
  ) %>%
  mutate(
    median_difference = LateProgressor_Score1 - LateNonProgressor_Score1
  ) %>%
  arrange(desc(abs(median_difference)))


effect_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(
    median_prog = median(score[signature=="LateProgressor_Score1"]),
    median_np = median(score[signature=="LateNonProgressor_Score1"]),
    median_diff = median_prog - median_np,
    mad_all = mad(score),
    standardized_effect = median_diff/mad_all
  )


## Combined
diff_violin_df <- NOD_T1D_Timepoints@meta.data %>%
  dplyr::select(
    celltype,
    EarlyProgressor_Score1,
    EarlyNonProgressor_Score1,
    IntermediateProgressor_Score1,
    IntermediateNonProgressor_Score1,
    LateProgressor_Score1,
    LateNonProgressor_Score1
  ) %>%
  filter(!is.na(celltype)) %>%
  mutate(
    Early = EarlyProgressor_Score1 - EarlyNonProgressor_Score1,
    Intermediate = IntermediateProgressor_Score1 - IntermediateNonProgressor_Score1,
    Late = LateProgressor_Score1 - LateNonProgressor_Score1
  ) %>%
  dplyr::select(celltype, Early, Intermediate, Late) %>%
  pivot_longer(
    cols = c(Early, Intermediate, Late),
    names_to = "stage",
    values_to = "score_difference"
  ) %>%
  mutate(
    stage = factor(stage, levels = c("Early", "Intermediate", "Late"))
  )

ggplot(diff_violin_df,
       aes(x = celltype,
           y = score_difference,
           fill = stage)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.75,
    position = position_dodge(width = 0.8)
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Progressor - Non-Progressor Score Difference",
    title = "P vs NP DEG Score Difference Across NOD Pancreas Cell Types",
    fill = NULL
  ) +
  scale_fill_manual(
    values = c(
      "Early" = "#F4A6A6",
      "Intermediate" = "#E60000",
      "Late" = "#8B0000"
    )
  ) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )

# Difference in module score

#Early
NOD_T1D_Timepoints$EarlyPvsNP <-
  NOD_T1D_Timepoints$EarlyProgressor_Score1 -
  NOD_T1D_Timepoints$EarlyNonProgressor_Score1

bar_df <- NOD_T1D_Timepoints@meta.data %>%
  dplyr::select(celltype, EarlyPvsNP) %>%
  filter(!is.na(celltype)) %>%
  group_by(celltype) %>%
  summarise(
    median_enrichment = median(EarlyPvsNP, na.rm = TRUE),
    mad_enrichment = mad(EarlyPvsNP, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_enrichment)) %>%
  mutate(
    celltype = factor(celltype, levels = celltype),
    direction = ifelse(median_enrichment >= 0,
                       "Progressor Score",
                       "Non-Progressor Score")
  )
  
ggplot(bar_df,
       aes(x = celltype,
           y = median_enrichment,
           fill = direction)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             linewidth = 0.7) +
  coord_flip() +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Median Progressor vs Non-Progressor Score Difference",
    title = "Early Stage DEGs-NOD Pancreas"
  ) +
  scale_fill_manual(
    values = c(
      "Progressor Score" = "#E60000",
      "Non-Progressor Score" = "#3D7A60"
    )
  ) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 18),
    legend.position = "none"
  )


#Intermediate
NOD_T1D_Timepoints$IntermediatePvsNP <-
  NOD_T1D_Timepoints$IntermediateProgressor_Score1 -
  NOD_T1D_Timepoints$IntermediateNonProgressor_Score1

bar_df <- NOD_T1D_Timepoints@meta.data %>%
  dplyr::select(celltype, IntermediatePvsNP) %>%
  filter(!is.na(celltype)) %>%
  group_by(celltype) %>%
  summarise(
    median_enrichment = median(IntermediatePvsNP, na.rm = TRUE),
    mad_enrichment = mad(IntermediatePvsNP, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_enrichment)) %>%
  mutate(
    celltype = factor(celltype, levels = celltype),
    direction = ifelse(median_enrichment >= 0,
                       "Progressor Score",
                       "Non-Progressor Score")
  )

ggplot(bar_df,
       aes(x = celltype,
           y = median_enrichment,
           fill = direction)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             linewidth = 0.7) +
  coord_flip() +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Median Progressor vs Non-Progressor Score Difference",
    title = "Intermediate Stage DEGs-NOD Pancreas"
  ) +
  scale_fill_manual(
    values = c(
      "Progressor Score" = "#E60000",
      "Non-Progressor Score" = "#3D7A60"
    )
  ) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 18),
    legend.position = "none"
  )

#Late
NOD_T1D_Timepoints$LatePvsNP <-
  NOD_T1D_Timepoints$LateProgressor_Score1 -
  NOD_T1D_Timepoints$LateNonProgressor_Score1

bar_df <- NOD_T1D_Timepoints@meta.data %>%
  dplyr::select(celltype, LatePvsNP) %>%
  filter(!is.na(celltype)) %>%
  group_by(celltype) %>%
  summarise(
    median_enrichment = median(LatePvsNP, na.rm = TRUE),
    mad_enrichment = mad(LatePvsNP, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_enrichment)) %>%
  mutate(
    celltype = factor(celltype, levels = celltype),
    direction = ifelse(median_enrichment >= 0,
                       "Progressor Score",
                       "Non-Progressor Score")
  )

ggplot(bar_df,
       aes(x = celltype,
           y = median_enrichment,
           fill = direction)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             linewidth = 0.7) +
  coord_flip() +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Median Progressor vs Non-Progressor Score Difference",
    title = "Late Stage DEGs-NOD Pancreas"
  ) +
  scale_fill_manual(
    values = c(
      "Progressor Score" = "#E60000",
      "Non-Progressor Score" = "#3D7A60"
    )
  ) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 18),
    legend.position = "none"
  )

#Combined

# Summarize median score difference by cell type and stage
bar_df_all <- NOD_T1D_Timepoints@meta.data %>%
  dplyr::select(celltype, EarlyPvsNP, IntermediatePvsNP, LatePvsNP) %>%
  filter(!is.na(celltype)) %>%
  pivot_longer(
    cols = c(EarlyPvsNP, IntermediatePvsNP, LatePvsNP),
    names_to = "stage",
    values_to = "score_difference"
  ) %>%
  mutate(
    stage = recode(
      stage,
      "EarlyPvsNP" = "Early",
      "IntermediatePvsNP" = "Intermediate",
      "LatePvsNP" = "Late"
    ),
    stage = factor(stage, levels = c("Early", "Intermediate", "Late"))
  ) %>%
  group_by(stage, celltype) %>%
  summarise(
    median_score_difference = median(score_difference, na.rm = TRUE),
    mad_score_difference = mad(score_difference, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    direction = ifelse(
      median_score_difference >= 0,
      "Progressor Score",
      "Non-Progressor Score"
    )
  )

ggplot(bar_df_all,
       aes(x = celltype,
           y = median_score_difference,
           fill = direction)) +
  geom_col(width = 0.75) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  coord_flip() +
  facet_wrap(~ stage, nrow = 1) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Median Progressor vs Non-Progressor Score Difference",
    title = "Progressor vs Non-Progressor DEG Scores Across NOD Pancreas Cell Types"
  ) +
  scale_fill_manual(
    values = c(
      "Progressor Score" = "#E60000",
      "Non-Progressor Score" = "#3D7A60"
    )
  ) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.position = "none"
  )

#4. Map Top 100 Genes ----
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

T1D_Blood = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/Blood/T1D_Seurat_Object_Final_v2.rds")

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

T1D_Spleen = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/MAI_T1D/spleen_result_annotated_v2.rds")
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


