
# Load the Libraries ----

pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, dnapath, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)

library(ggplot2)
library(dplyr)
devtools::install_github("mixOmicsTeam/mixOmics")
library(mixOmics)
packageVersion("mixOmics")  # Should be 6.24.0 or similar
# Then unload conflicting packages if loaded
detach("package:pls", unload = TRUE)
detach("package:caret", unload = TRUE)
detach("package:plsRglm", unload = TRUE)
library(colorRamp2)
library(ComplexHeatmap)
install.packages('grr')
library(grr)
# Load libraries
BiocManager::install("scater")
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)

devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.7.tar.gz")
library(Matrix.utils)##
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)


# Import Data ----
#Metadata Importing
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Pathway Analysis/metadata_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Pathway Analysis/metadata_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Pathway Analysis/metadata_MetabolomicsCohort_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)

# Merge metadata by columns 
meta_combined <- rbind(meta_batch1, meta_batch2,meta_batch3)
# Preview the combined metadata
head(meta_combined)

#Counts Data Importing
counts_batch1 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Pathway Analysis/gene_expected_count.annot_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch1 <- na.omit(counts_batch1)

counts_batch2 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Pathway Analysis/gene_expected_count.annot_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch2 <- na.omit(counts_batch2)

counts_batch3 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Pathway Analysis/gene_expected_count.annot_MetabolomicsCohort_Week6.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch3 <- na.omit(counts_batch3)


#Remove duplicate names
counts_batch1 <- counts_batch1[!duplicated(counts_batch1[, 1]), ]
genes <- counts_batch1[, 1]
rownames(counts_batch1) <- genes
counts_batch1 <- counts_batch1[, -1]

counts_batch2 <- counts_batch2[!duplicated(counts_batch2[, 1]), ]
genes <- counts_batch2[, 1]
rownames(counts_batch2) <- genes
counts_batch2 <- counts_batch2[, -1]

counts_batch3 <- counts_batch3[!duplicated(counts_batch3[, 1]), ]
genes <- counts_batch3[, 1]
rownames(counts_batch3) <- genes
counts_batch3 <- counts_batch3[, -1]

#Combine data
# Merge counts_batch1 and counts_batch2
combined_counts <- merge(counts_batch1, counts_batch2, by = "row.names", all = TRUE)
# Merge the result with counts_batch3
combined_counts <- merge(combined_counts, counts_batch3, by.x = "Row.names", by.y = "row.names", all = TRUE)

# Set rownames back to genes
rownames(combined_counts) <- combined_counts$Row.names
combined_counts <- combined_counts[, -1]

# Preview the combined dataset
head(combined_counts)

# Batch Corrct Data----

early_data <- combined_counts[, meta_combined$Time == "Early"]
early_data <- na.omit(early_data)
meta_early <- meta_combined[meta_combined$Time == "Early", ]


early_data_Batchcorrected <- flexiDEG.function1(early_data, meta_early, # Run Function 1
                                                convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                                batches = T, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
early_data_Batchcorrected <- early_data_Batchcorrected[!grepl("^Gm[0-9]", rownames(early_data_Batchcorrected)), ]

#Import DESeq Results for Early Stage
Early_DEResults<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Early_NOD_ProgressorVsNonProgressor.csv")
rownames(Early_DEResults)<-Early_DEResults$X


# PLSDA-Early Stage DEGs----
#Subset Significant genes for Early Stag
Early_DEG_genes <- rownames(Early_DEResults)[ abs(Early_DEResults$log2FoldChange) >=1 &
  Early_DEResults$pvalue <= 0.05 
]

# Subset the sig genes from the batch-corrected data
early_data_Batchcorrected_DEGFiltered <- early_data_Batchcorrected[rownames(early_data_Batchcorrected) %in% Early_DEG_genes, ]



set.seed(123)  # For reproducibility
# Prepare data
X <- t(early_data_Batchcorrected_DEGFiltered)  # samples x genes
Y <- meta_early$Group                          # class labels
Y <- as.factor(meta_early$Group)  # Convert class labels to factor

# Run PLS-DA
plsda_model <-  mixOmics::plsda(X, Y, ncomp = 2)

# Extract scores
scores_top <- plsda_model$variates$X
# Plot PLSDA with Top 100 Genes
plot_df_top <- data.frame(LV1 = scores_top[,1],
                          LV2 = scores_top[,2],
                          Group = Y)
# Define colors & shapes
group_colors <- c("Progressor" = "#E60000", "Non-Progressor" = "#3D7A60")
group_shapes <- c("Progressor" = 16, "Non-Progressor" = 17)

ggplot(plot_df_top, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.68) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "PLS-DA: Early Stage DEG Genes ",
       x = "Latent Variable 1", y = "Latent Variable 2") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )




# PLSDA-Top 100 Genes----

#Subset Significant genes for Early Stag
sig_genes <- rownames(Early_DEResults)[
  Early_DEResults$pvalue <= 0.05 
]

# Subset the sig genes from the batch-corrected data
early_data_Batchcorrected_Filtered <- early_data_Batchcorrected[rownames(early_data_Batchcorrected) %in% sig_genes, ]


set.seed(123)  # For reproducibility
# Prepare data
X <- t(early_data_Batchcorrected_Filtered)  # samples x genes
Y <- meta_early$Group                          # class labels
Y <- as.factor(meta_early$Group)  # Convert class labels to factor

# Run PLS-DA
plsda_model <-  mixOmics::plsda(X, Y, ncomp = 2)

# Extract VIP scores for all components
vip_scores <- vip(plsda_model)
# Extract VIPs from component 1
vip_lv1 <- vip_scores[, 1]
# Export VIP scores to CSV
write.csv(vip_scores, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Prediction Analysis/PLSDA_Vip_scores_new.csv", row.names = TRUE)

top_gene_names <- names(sort(vip_scores[,1], decreasing = TRUE)[1:100])
# Subset X to top VIP genes only
X_top <- X[, top_gene_names]
# Rerun PLS-DA using top 100 Genes (For Plotting)
plsda_top_model <- mixOmics::plsda(X_top, Y, ncomp = 2)
# Extract scores
scores_top <- plsda_top_model$variates$X
# Plot PLSDA with Top 100 Genes
plot_df_top <- data.frame(LV1 = scores_top[,1],
                          LV2 = scores_top[,2],
                          Group = Y)
# Define colors & shapes
group_colors <- c("Progressor" = "#E60000", "Non-Progressor" = "#3D7A60")
group_shapes <- c("Progressor" = 16, "Non-Progressor" = 17)

ggplot(plot_df_top, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.68) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "PLS-DA: Top 100 VIP Genes ",
       x = "Latent Variable 1", y = "Latent Variable 2") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )

# Heatmap-Top 100 Genes----
early_data_Batchcorrected_Top100<- early_data_Batchcorrected[rownames(early_data_Batchcorrected) %in% top_gene_names, ]
# Expression matrix (genes x samples)
mat <- as.matrix(early_data_Batchcorrected_Top100)

# Metadata: group info for samples
group <- meta_early$Group
names(group) <- meta_early$Samples

# Annotation colors
group_col <- c("Progressor" = "#E60000", "Non-Progressor" = "#3D7A60")

# Column annotation
ha_col <- HeatmapAnnotation(
  Group = group[colnames(mat)],
  col = list(Group = group_col),
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  annotation_legend_param = list(title = "Group")
)

# Row-wise Z-score scaling
mat_scaled <- t(scale(t(mat)))

# Color function for heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))


# 1) choose which genes you want to label
genes_to_label <- c(
  # Immune / inflammation / T cell regulation
  "Pdcd1", "Cd40lg", "Ccr8", "Cpvl", "Cfi", "Tnfaip8l1",
  #Metabolic
  "Cpt2", "Gpd2", "Aacs", "Mpst",
  #Chromatin Regulators
  "Ncor1", "Sp8","Yy2","Zfp703",
  #Trafficking Proteins
  "Nedd4l","Psme3", "Rab3c","Snx3",
  #Stromal
  "Col10a1","Bambi","Creld1",
  #ERV
  "Zcchc18",
  #Othrs
  "Tm9sf4", "Cdkn2c", "Stmn1","Eogt","Lbh","Rbp4",
  "Dusp12", "Prkcz2", "Gsk3a", "Wnk1",  
  "Pemt",  "Cth",  
  "Ephx3", "Cnnm2", "Pgrmc2","Cdc26","Nppa"
)
at_idx <- which(rownames(mat_scaled) %in% genes_to_label)
lab    <- rownames(mat_scaled)[at_idx]

# 2) build the row "mark" annotation (leader lines + labels)
ra <- rowAnnotation(
  mark = anno_mark(
    at = at_idx,
    labels = lab,
    side = "right",
    labels_gp = gpar(fontsize = 12, fontface = "italic"),
    link_width = unit(10, "mm"),
    link_height = unit(2, "mm")
  )
)

# Heatmap with gene names
ht <- Heatmap(
  mat_scaled,
  name = "Z-score",
  top_annotation = ha_col,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 12, fontface = "italic"),  # Italicized gene names
  show_column_names = FALSE,
  heatmap_legend_param = list(title = "Expression", legend_direction = "vertical")
)

pdf("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Figure/T1D_Top100Genes_heatmap.pdf", width = 7, height = 8)
draw(ht + ra, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


write.csv(mat_scaled,"/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Result/T1D_Top100Genes_Heatmap_scaled.csv")

out_file <- "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Result/T1D_Top100Genes_Heatmap_scaled.csv"

# group vector aligned to mat_scaled columns
group_row <- group[colnames(mat_scaled)]

con <- file(out_file, open = "wt")

# 1) header row: blank cell for gene names + sample names
writeLines(
  paste(c("", colnames(mat_scaled)), collapse = ","),
  con
)
# 2) second row: label + group values
writeLines(
  paste(c("Group", as.character(group_row)), collapse = ","),
  con
)
# 3) data (genes x samples), with gene names in first column
write.table(
  mat_scaled,
  file = con,
  sep = ",",
  quote = FALSE,
  col.names = FALSE,
  row.names = TRUE
)
close(con)


# Heatmap-Early DEG Genes----

# Expression matrix (genes x samples)
mat <- as.matrix(early_data_Batchcorrected_DEGFiltered)

# Metadata: group info for samples
group <- meta_early$Group
names(group) <- meta_early$Samples

# Annotation colors
group_col <- c("Progressor" = "#E60000", "Non-Progressor" = "#3D7A60")

# Column annotation
ha_col <- HeatmapAnnotation(
  Group = group[colnames(mat)],
  col = list(Group = group_col),
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  annotation_legend_param = list(title = "Group")
)

# Row-wise Z-score scaling
mat_scaled <- t(scale(t(mat)))
library(ComplexHeatmap)
library(circlize)

# Color function for heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))




# Heatmap with gene names
 Heatmap(
  mat_scaled,
  name = "Z-score",
  top_annotation = ha_col,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),  # Italicized gene names
  show_column_names = FALSE,
  heatmap_legend_param = list(title = "Expression", legend_direction = "vertical")
)



write.csv(mat_scaled,"/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Result/T1D_Top100Genes_Heatmap_scaled.csv")

out_file <- "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Result/T1D_Top100Genes_Heatmap_scaled.csv"

# group vector aligned to mat_scaled columns
group_row <- group[colnames(mat_scaled)]

con <- file(out_file, open = "wt")

# 1) header row: blank cell for gene names + sample names
writeLines(
  paste(c("", colnames(mat_scaled)), collapse = ","),
  con
)
# 2) second row: label + group values
writeLines(
  paste(c("Group", as.character(group_row)), collapse = ","),
  con
)
# 3) data (genes x samples), with gene names in first column
write.table(
  mat_scaled,
  file = con,
  sep = ",",
  quote = FALSE,
  col.names = FALSE,
  row.names = TRUE
)
close(con)

# Top 100 Genes-Human ----
T1DGene_top100_humanorthologs<- read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Input/SVC Model/Top100Genes_T1D_HumanOrthologs.csv", sep=",", header=T,check.names = FALSE)
Top100_Human<- T1DGene_top100_humanorthologs$human_symbol
Top100_Human <- Top100_Human[!(Top100_Human=='N/A')]


## Human Spleen----

T1D_Spleen = readRDS("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/SingleCellRNASeq/MAI_T1D/spleen_result_annotated_v2.rds")
T1D_Spleen <- RenameIdents(
  T1D_Spleen,
  "B" = "B Cell",
  "T" = "T Cell",
  "Cycling" = "Cycling Cell",
  "NK"="NK Cell"
)

unique(Idents(T1D_Spleen))
unique(T1D_Spleen$group)

###  Psuedobuk Mean expression (Across CellTypes_-----

DefaultAssay(T1D_Spleen) <- "RNA"
counts_spleen <- GetAssayData(T1D_Spleen, assay = "RNA", layer = "counts")  # genes x cells
sample_id_spleen <- factor(T1D_Spleen$samples)

# 1) Sum counts per sample (samples x genes)
pb_sample_sum_spleen <- aggregate.Matrix(
  t(counts_spleen),
  groupings = data.frame(sample_id = sample_id_spleen),
  fun = "sum"
)

# 2) Convert to mean per cell in sample = (sum counts) / (number of cells in sample)
n_cells_per_sample_spleen <- as.numeric(table(sample_id_spleen))  # aligned to levels(sample_id_spleen)

pb_sample_mean_spleen <- pb_sample_sum_spleen
pb_sample_mean_spleen[] <- pb_sample_mean_spleen / n_cells_per_sample_spleen   # row-wise divide

# Convert to genes x samples to match your previous object naming
bulk_weighted_spleen <- t(pb_sample_mean_spleen)   # genes x samples
colnames(bulk_weighted_spleen) <- levels(sample_id_spleen)


### Plot PLSDA Uisng Top 100 Genes----
genes_use_spleen <- intersect(Top100_Human, rownames(bulk_weighted_spleen))
length(genes_use_spleen) #93 Genes
bulk_top100_spleen <- bulk_weighted_spleen[genes_use_spleen, ,drop = FALSE]

X_spleen <- t(bulk_top100_spleen)   # samples x genes
meta_spleen <- data.frame(
  Sample = colnames(bulk_weighted_spleen),
  Group  = ifelse(grepl("__T1D", colnames(bulk_weighted_spleen)), "T1D", "ND")
)

Y_spleen <- factor(meta_spleen$Group)
set.seed(123)

plsda_model_spleen <- mixOmics::plsda(X_spleen, Y_spleen, ncomp = 2,scale = TRUE)
scores_spleen <- plsda_model_spleen$variates$X

plot_df_spleen <- data.frame(
  LV1 = scores_spleen[,1],
  LV2 = scores_spleen[,2],
  Group = Y_spleen
)

group_colors <- c("T1D" = "#E60000", "ND" = "#3D7A60")
group_shapes <- c("T1D" = 16, "ND" = 17)


ggplot(plot_df_spleen, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.68) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "Human Spleen",
       x = "Latent Variable 1", y = "Latent Variable 2") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )



# Human pLN----

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

unique(Idents(T1D_pln))
unique(T1D_pln$group)

###  Psuedobuk Mean expression (Across CellTypes_-----

DefaultAssay(T1D_pln) <- "RNA"
T1D_pln[["RNA"]] <- JoinLayers(T1D_pln[["RNA"]])
counts_pln <- GetAssayData(T1D_pln, assay = "RNA", layer = "counts")  # genes x cells
sample_id_pln <- factor(T1D_pln$samples)

# 1) Sum counts per sample (samples x genes)
pb_sample_sum_pln <- aggregate.Matrix(
  t(counts_pln),
  groupings = data.frame(sample_id = sample_id_pln),
  fun = "sum"
)

# 2) Convert to mean per cell in sample = (sum counts) / (number of cells in sample)
n_cells_per_sample_pln <- as.numeric(table(sample_id_pln))  # aligned to levels(sample_id_pln)

pb_sample_mean_pln <- pb_sample_sum_pln
pb_sample_mean_pln[] <- pb_sample_sum_pln / n_cells_per_sample_pln   # row-wise divide

# Convert to genes x samples to match your previous object naming
bulk_weighted_pln <- t(pb_sample_mean_pln)   # genes x samples
colnames(bulk_weighted_pln) <- levels(sample_id_pln)


### Plot PLSDA Uisng Top 100 Genes----
genes_use_pln <- intersect(Top100_Human, rownames(bulk_weighted_pln))
length(genes_use_pln) #93
bulk_top100_pln <- bulk_weighted_pln[genes_use_pln, ,drop=FALSE]

X_pln <- t(bulk_top100_pln)   # samples x genes
meta_pln <- data.frame(
  Sample = colnames(bulk_weighted_pln),
  Group  = ifelse(grepl("__T1D", colnames(bulk_weighted_pln)), "T1D", "ND")
)

Y_pln <- factor(meta_pln$Group)
set.seed(123)

plsda_model_pln <- mixOmics::plsda(X_pln, Y_pln, ncomp = 2,scale = TRUE)
scores_pln <- plsda_model_pln$variates$X

plot_df_pln <- data.frame(
  LV1 = scores_pln[,1],
  LV2 = scores_pln[,2],
  Group = Y_pln
)

group_colors <- c("T1D" = "#E60000", "ND" = "#3D7A60")
group_shapes <- c("T1D" = 16, "ND" = 17)

ggplot(plot_df_pln, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.68) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "Human pLN",
       x = "Latent Variable 1", y = "Latent Variable 2") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )


# Human Blood----


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


###  Psuedobuk Mean expression (Across CellTypes_-----

DefaultAssay(T1D_Blood) <- "RNA"
counts_Blood <- GetAssayData(T1D_Blood, assay = "RNA", layer = "counts")  # genes x cells
sample_id_Blood <- factor(T1D_Blood$Sample_ID)

# 1) Sum counts per sample (samples x genes)
pb_sample_sum_Blood <- aggregate.Matrix(
  t(counts_Blood),
  groupings = data.frame(sample_id = sample_id_Blood),
  fun = "sum"
)

# 2) Convert to mean per cell in sample = (sum counts) / (number of cells in sample)
n_cells_per_sample_Blood <- as.numeric(table(sample_id_Blood))  # aligned to levels(sample_id_Blood)

pb_sample_mean_Blood <- pb_sample_sum_Blood
pb_sample_mean_Blood[] <- pb_sample_sum_Blood / n_cells_per_sample_Blood   # row-wise divide

# Convert to genes x samples to match your previous object naming
bulk_weighted_Blood <- t(pb_sample_mean_Blood)   # genes x samples
colnames(bulk_weighted_Blood) <- levels(sample_id_Blood)


### Plot PLSDA Uisng Top 100 Genes----
genes_use_Blood <- intersect(Top100_Human, rownames(bulk_weighted_Blood))
length(genes_use_Blood)
bulk_top100_Blood <- bulk_weighted_Blood[genes_use_Blood, ,drop=FALSE]
#73 Genes
X_Blood <- t(bulk_top100_Blood)   # samples x genes
meta_Blood <- data.frame(
  Sample = colnames(bulk_weighted_Blood),
  Group  = ifelse(grepl("T1D-", colnames(bulk_weighted_Blood)), "T1D", "ND")
)

Y_Blood <- factor(meta_Blood$Group)
set.seed(123)

plsda_model_Blood <- mixOmics::plsda(X_Blood, Y_Blood, ncomp = 2,scale = TRUE)
scores_Blood <- plsda_model_Blood$variates$X

plot_df_Blood <- data.frame(
  LV1 = scores_Blood[,1],
  LV2 = scores_Blood[,2],
  Group = Y_Blood
)

group_colors <- c("T1D" = "#E60000", "ND" = "#3D7A60")
group_shapes <- c("T1D" = 16, "ND" = 17)

ggplot(plot_df_Blood, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.68) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "Human Blood",
       x = "Latent Variable 1", y = "Latent Variable 2") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  )

