pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, dnapath, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library(dplyr)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(VennDiagram)
library(grid)



# --- Final combined TERM2GENE data frame ---
mm_all_df <- rbind(mm_h_df, mm_kegg_df)

getwd()

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

# Early Stage DESEQ----

early_data <- combined_counts[, meta_combined$Time == "Early"]
early_data <- na.omit(early_data)
meta_early <- meta_combined[meta_combined$Time == "Early", ]

early_data <- flexiDEG.function1(early_data, meta_early, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
early_data <- early_data[!grepl("^Gm[0-9]", rownames(early_data)), ]


# DESeq2 analysis for Early Time
dds_early <- DESeqDataSetFromMatrix(countData = early_data, colData = meta_early, design = ~ Batch+Group)
dds_early <- DESeq(dds_early)
results_early <- as.data.frame(results(dds_early, contrast = c("Group", "Progressor", "Non-Progressor")))

# Replace NA p values with 1
results_early$pvalue[is.na(results_early$pvalue)] <- 1
results_early$padj[is.na(results_early$padj)] <- 1



# Define color mapping with stronger contrast for visibility
keyvals <- ifelse(results_early$log2FoldChange > 1 & results_early$pvalue < 0.05, "#FF0000",  # Bright red for nominally significant up
                                ifelse(results_early$log2FoldChange < -1 & results_early$pvalue < 0.05, "#1E90FF",  # Bright blue for nominally significant down
                                       "gray"))  # Gray for non-significant genes

# Assign names for legend
names(keyvals) <- ifelse(results_early$log2FoldChange > 1 & results_early$pvalue < 0.05, "p-value < 0.05 & log[2]FC > 1",
                                       ifelse(results_early$log2FoldChange < -1 & results_early$pvalue < 0.05, "p-value < 0.05 & log[2]FC < -1",
                                              "Not Significant"))

point_sizes <-  2.0

# Genes to label
genes_to_label_early <- c(
  "Cxcl15","Cfd","Retn",
  "Cd40lg","Cd6","Ctla4","Ms4a1","Gata1",
  "Ccr8",       
  "Fgg",
  "Ly6G",
  "Ctla4",
  "Cd163l1",   
  "Cpvl",      
  "Il23r",      
  "Lgals12",    
  "Retnlg",     
  "S100a8",     
  "S100a9",     
  "Vnn3"       
)





EnhancedVolcano(results_early,
                lab = rownames(results_early),
                selectLab = genes_to_label_early,
                x = 'log2FoldChange',
                y = 'pvalue',  
                title = 'Early Timepoint',
                pCutoff = 0.05,  
                FCcutoff = 1, 
                pointSize = point_sizes,  
                labSize = 8,
                colAlpha = 0.75,
                legendLabels = c("Not Significant", "p-value < 0.05", "FDR < 0.1", "Both"),
                legendPosition = 'right',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                boxedLabels = TRUE,
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"), 
                ylab = expression("-log"[10] ~ "(p-value)"),  
                xlim = c(-8, 8),  # Set x-axis from -8 to 8
                ylim = c(0, 6),  # Set y-axis limit to 6
                max.overlaps = 40,
                colCustom = keyvals, 
                caption = NULL  
)



# Save the data frame as CSV
write.csv(results_early, file = "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Early_NOD_ProgressorVsNonProgressor.csv")


#Intermediate Stage DESEQ----

Intermediate_data <- combined_counts[, meta_combined$Time == "Intermediate"]
Intermediate_data <- na.omit(Intermediate_data)
meta_Intermediate <- meta_combined[meta_combined$Time == "Intermediate", ]

Intermediate_data <- flexiDEG.function1(Intermediate_data, meta_Intermediate, # Run Function 1
                                        convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                        batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
Intermediate_data <- Intermediate_data[!grepl("^Gm[0-9]", rownames(Intermediate_data)), ]


# DESeq2 analysis for Intermediate Time
dds_Intermediate <- DESeqDataSetFromMatrix(countData = Intermediate_data, colData = meta_Intermediate, design = ~ Batch+Group)
dds_Intermediate <- DESeq(dds_Intermediate)
results_Intermediate <- as.data.frame(results(dds_Intermediate, contrast = c("Group", "Progressor", "Non-Progressor")))

# Replace NA p-values with 1
results_Intermediate$pvalue[is.na(results_Intermediate$pvalue)] <- 1
results_Intermediate$padj[is.na(results_Intermediate$padj)] <- 1

# Define color mapping with stronger contrast for visibility
keyvals <- ifelse(results_Intermediate$log2FoldChange > 1 & results_Intermediate$pvalue < 0.05, "#FF0000",  # Bright red for nominally significant up
                                ifelse(results_Intermediate$log2FoldChange < -1 & results_Intermediate$pvalue < 0.05, "#1E90FF",  # Bright blue for nominally significant down
                                       "gray"))  # Gray for non-significant genes

# Assign names for legend
names(keyvals) <- ifelse(results_Intermediate$log2FoldChange > 1 & results_Intermediate$pvalue < 0.05, "p-value < 0.05 & log[2]FC > 1",
                                       ifelse(results_Intermediate$log2FoldChange < -1 & results_Intermediate$pvalue < 0.05, "p-value < 0.05 & log[2]FC < -1",
                                              "Not Significant"))
point_sizes <-  2.0

# Genes to label

genes_to_label_intermediate <- c("Cd300ld5", "Mlf1", "Mmp3",
                                 "Fpr1",
                                 "Retn",     
                                 "Defb1",   
                                 "Ifi44",   
                                 "Cxcl15",   
                                 "Rtp4",     
                                 "Ido1",     
                                 "Gzmk",     
                                 "Xcl1",      
                                 "Cd4",    
                                 "Cd7",     
                                 "Cd69",     
                                 "Rorc",    
                                 "Ighg1",    
                                 "Xcl1",     
                                 "Fcrl6",   
                                 "Fpr1",    
                                 "Ifit1",    
                                 "Ifit3",     
                                 "Klra4",  
                                 "Klra5", 
                                 "Ifi208",
                                 "Btnl9",
                                 "Ido1",
                                 "Rorc"
)


EnhancedVolcano(results_Intermediate,
                lab = rownames(results_Intermediate),
                selectLab = genes_to_label_intermediate,
                x = 'log2FoldChange',
                y = 'pvalue',  
                title = 'Intermediate Timepoint',
                pCutoff = 0.05,  
                FCcutoff = 1, 
                pointSize = point_sizes, 
                labSize = 8,
                colAlpha = 0.75,
                legendLabels = c("Not Significant", "p-value < 0.05", "FDR < 0.1", "Both"),
                legendPosition = 'right',
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"),  
                ylab = expression("-log"[10] ~ "(p-value)"),  
                xlim = c(-8, 8),  
                ylim = c(0, 6),  
                colCustom = keyvals,  
                caption = NULL 
)


# Save the data frame as CSV
write.csv(results_Intermediate, file = "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Intermediate_NOD_ProgressorVsNonProgressor.csv")

# Late Stage DESEQ----

Late_data <- combined_counts[, meta_combined$Time == "Late"]
Late_data <- na.omit(Late_data)
meta_Late <- meta_combined[meta_combined$Time == "Late", ]

Late_data <- flexiDEG.function1(Late_data, meta_Late, # Run Function 1
                                convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
Late_data <- Late_data[!grepl("^Gm[0-9]", rownames(Late_data)), ]


# DESeq2 analysis for Late Time
dds_Late <- DESeqDataSetFromMatrix(countData = Late_data, colData = meta_Late, design = ~ Batch+Group)
dds_Late <- DESeq(dds_Late)
results_Late <- as.data.frame(results(dds_Late, contrast = c("Group", "Progressor", "Non-Progressor")))
# Replace NA p values with 1
results_Late$padj[is.na(results_Late$padj)] <- 1
results_Late$pvalue[is.na(results_Late$pvalue)] <- 1


# Define color mapping with enhanced contrast
keyvals <- ifelse(results_Late$log2FoldChange > 1 & results_Late$pvalue < 0.05, "#FF0000",  # Bright red for nominally significant up
                                ifelse(results_Late$log2FoldChange < -1 & results_Late$pvalue < 0.05, "#1E90FF",  # Bright blue for nominally significant down
                                       "gray"))  # Gray for non-significant genes

# Assign names for legend
names(keyvals) <- ifelse(results_Late$log2FoldChange > 1 & results_Late$pvalue < 0.05, "p-value < 0.05 & log[2]FC > 1",
                                       ifelse(results_Late$log2FoldChange < -1 & results_Late$pvalue < 0.05, "p-value < 0.05 & log[2]FC < -1",
                                              "Not Significant"))

point_sizes <-  2.0

# Genes to label
genes_to_label_late <- c("Ighg2c", "Il31ra", "Tafa2","Cd7", "Il13", "Il17b", "Il17re", "Prf1", "Rorc", "Trbc2", "Tnfrsf17", "Trcg1", "Cd163",    # Classic M2 macrophage scavenger receptor
                         "Cd209b",   
                         "Vnn1",
                         "Rorc",
                         "Cxcl17", 
                         "Ccl17",
                         "Ccr8",
                         "Cxcl3",
                         "Tox" 
)  



EnhancedVolcano(results_Late,
                lab = rownames(results_Late),
                selectLab = genes_to_label_late,
                x = 'log2FoldChange',
                y = 'pvalue',  
                title = 'Late Timepoint',
                boxedLabels = TRUE,
                pCutoff = 0.05,  
                FCcutoff = 1,  
                pointSize = point_sizes,  
                labSize = 8,
                colAlpha = 0.75,
                legendLabels = c("Not Significant", "p-value < 0.05", "FDR < 0.1", "Both"),
                legendPosition = 'right',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps=30,
                xlab = expression("log"[2] ~ "Fold Change (Progressor / Non-Progressor)"), 
                ylab = expression("-log"[10] ~ "(p-value)"), 
                xlim = c(-8, 8),  
                ylim = c(0, 6), 
                colCustom = keyvals,  
                caption = NULL  
)



# Save the data frame as CSV
write.csv(results_Late, file = "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Late_NOD_ProgressorVsNonProgressor.csv")

#GSEA Analysis----

msigdbr_collections() 

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

#### Early----

results_early<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Early_NOD_ProgressorVsNonProgressor.csv",row.names = 1)

#Build ranked gene lists using DESeq2 Wald stat
lfc_vector_early  <- results_early$stat;  names(lfc_vector_early)  <- rownames(results_early)

# Drop NAs
lfc_vector_early  <- lfc_vector_early[!is.na(lfc_vector_early)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_early  <- sort(lfc_vector_early,  decreasing = TRUE)


gsea_results_early <- GSEA(
  geneList = lfc_vector_early, 
  minGSSize = 5, 
  maxGSSize = 500, 
  pvalueCutoff = 1,
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = mm_all_df  
)

# Extract results for Early
gsea_results_early_df <- as.data.frame(gsea_results_early)


#### Intermediate----
results_Intermediate<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Intermediate_NOD_ProgressorVsNonProgressor.csv",row.names = 1)

#Build ranked gene lists using DESeq2 Wald stat
lfc_vector_intermediate  <- results_Intermediate$stat;  names(lfc_vector_intermediate)  <- rownames(results_Intermediate)

# Drop NAs
lfc_vector_intermediate  <- lfc_vector_intermediate[!is.na(lfc_vector_intermediate)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_intermediate  <- sort(lfc_vector_intermediate,  decreasing = TRUE)


# Perform GSEA for Intermediate
gsea_results_intermediate <- GSEA(
  geneList = lfc_vector_intermediate, 
  minGSSize = 5, 
  maxGSSize = 500,
  pvalueCutoff = 1,
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = mm_all_df 
)

# Extract results for Intermediate
gsea_results_intermediate_df <- as.data.frame(gsea_results_intermediate)

#### Late -----
results_Late<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Late_NOD_ProgressorVsNonProgressor.csv",row.names = 1)

#Build ranked gene lists using DESeq2 Wald stat
lfc_vector_Late  <- results_Late$stat;  names(lfc_vector_Late)  <- rownames(results_Late)

# Drop NAs
lfc_vector_Late  <- lfc_vector_Late[!is.na(lfc_vector_Late)]

# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_Late  <- sort(lfc_vector_Late,  decreasing = TRUE)

# Perform GSEA for Late
gsea_results_late <- GSEA(
  geneList = lfc_vector_Late, 
  minGSSize = 5, 
  maxGSSize = 500, 
  pvalueCutoff = 1,
  eps = 0, 
  seed = TRUE,
  pAdjustMethod = "BH", 
  TERM2GENE = mm_all_df  
)

# Extract results for Late
gsea_results_late_df <- as.data.frame(gsea_results_late)

#### Combined GSEA----
#Add prefixes to column names to identify the source of each data frame
names(gsea_results_early_df)[-1] <- paste0("early_", names(gsea_results_early_df)[-1])
names(gsea_results_intermediate_df)[-1] <- paste0("intermediate_", names(gsea_results_intermediate_df)[-1])
names(gsea_results_late_df)[-1] <- paste0("late_", names(gsea_results_late_df)[-1])

# Merge the three data frames based on the 'Description' column
combined_gsea_df <- Reduce(function(x, y) merge(x, y, by = "ID"), 
                           list(gsea_results_early_df, gsea_results_intermediate_df, gsea_results_late_df))

# View the combined data frame
head(combined_gsea_df)
# Save the combined data frame to a CSV file
write.csv(combined_gsea_df, "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/GSEA_NODProgressorVsNonProgressor_AllTimepoints.csv", row.names = FALSE)

# Filter pathways where NES >= 1 or <= -1 AND pvalue < 0.05 for any stage
filtered_gsea_df <- combined_gsea_df[
  ((combined_gsea_df$early_NES >= 1 | combined_gsea_df$early_NES <= -1) & combined_gsea_df$early_p.adjust < 0.1) |
    ((combined_gsea_df$intermediate_NES >= 1 | combined_gsea_df$intermediate_NES <= -1) & combined_gsea_df$intermediate_p.adjust < 0.1) |
    ((combined_gsea_df$late_NES >= 1 | combined_gsea_df$late_NES <= -1) & combined_gsea_df$late_p.adjust < 0.1),
]
dim(filtered_gsea_df)
# View the filtered data
head(filtered_gsea_df)

# Save the filtered data frame to a CSV file
write.csv(filtered_gsea_df, "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/GSEA_Filtered_NODProgressorVsNonProgressor_AllTimepoints.csv", row.names = FALSE)



Immune_Metabolic_pathways <- c(
  # T cell-centric & adaptive immunity
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  # Myeloid / innate sensing & effector
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
  "KEGG_LYSOSOME",
  # Inflammation & complement
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "HALLMARK_APOPTOSIS",
  # Immune–metabolic crosstalk (β-cell stress, T cell metabolism)
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_INSULIN_SIGNALING_PATHWAY",
  "KEGG_PPAR_SIGNALING_PATHWAY",
  "KEGG_TRYPTOPHAN_METABOLISM",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
)

filtered_gsea_df<-as.data.frame(read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/GSEA_Filtered_NODProgressorVsNonProgressor_AllTimepoints.csv"))
# Reshape the data into long format for plotting
filtered_gsea_long <- filtered_gsea_df %>%
  select(ID, starts_with("early_"), starts_with("intermediate_"), starts_with("late_")) %>%
  pivot_longer(
    cols = -ID,
    names_to = c("timepoint", ".value"),
    names_pattern = "(early|intermediate|late)_(.+)"
  )

# Calculate bubble size as -log10(pvalue)
filtered_gsea_long <- filtered_gsea_long %>%
  mutate(bubble_size = -log10(p.adjust))

filtered_gsea_immune_top <- filtered_gsea_long %>%
  filter(ID %in% Immune_Metabolic_pathways)

filtered_gsea_immune_top <- filtered_gsea_long %>%
  filter(ID %in% Immune_Metabolic_pathways) %>%
  mutate(
    timepoint = case_when(
      timepoint == "early" ~ "Early",
      timepoint == "intermediate" ~ "Intermediate",
      timepoint == "late" ~ "Late",
      TRUE ~ timepoint
    ),
    ID = factor(ID, levels = Immune_Metabolic_pathways)  # arrange order here
  )


plot <- ggplot(filtered_gsea_immune_top, aes(x = timepoint, y = ID)) +
  geom_point(aes(size = bubble_size, color = NES)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(8, 30), name = bquote(bold("-log10(" * p[adj] * ")"))) +
  theme_minimal()+
  labs(
    x = "Time",
    y = "Pathway",
    color = bquote(bold("NES"))  # Bold NES 
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 40, face = "bold"),
    axis.text.y = element_text(size = 40),  
    axis.title.x = element_text(size = 60),  
    axis.title.y = element_text(size = 60),  
    legend.text = element_text(size = 30),  
    legend.title = element_text(size = 40, face = "bold"),  
    legend.position = "right",
    plot.title = element_text(size = 30, face = "bold")  
  )

plot


#Pan-Disease Analysis (EAE,Cancer,T1D)----

## DESeq----

##### T1D----
results_early<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Early_NOD_ProgressorVsNonProgressor.csv",row.names = 1)
results_early[is.na(results_early)] <- 1
# Filter for significant genes with log2FC >= 1 or <= -1 and p-value <= 0.05
T1D_DEGs_Early <- results_early[
  (results_early$log2FoldChange >= 1 | results_early$log2FoldChange <= -1) &
    results_early$pvalue <= 0.05,
]

##### Cancer 4T1----


#Metadata
meta_cancer <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/PanDiseaseInflammation/Cancer4T1_Metadata.csv", sep=",", header=T) # Metadata file
meta_cancer <- as.data.frame(meta_cancer)

#Counts Data
counts_cancer <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/PanDiseaseInflammation/Cancer4T1_RNASeq.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_cancer <- na.omit(counts_cancer)
counts_cancer <- counts_cancer[!duplicated(counts_cancer[, 1]), ]
genes <- counts_cancer[, 1]
rownames(counts_cancer) <- genes
counts_cancer <- counts_cancer[, -1]

counts_cancer <- flexiDEG.function1(counts_cancer, meta_cancer, # Run Function 1
                                    convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                    batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
counts_cancer <- counts_cancer[!grepl("^Gm[0-9]", rownames(counts_cancer)), ]

# DESEQ2
dds_cancer <- DESeqDataSetFromMatrix(countData = counts_cancer, colData = meta_cancer, design = ~ Group)
dds_cancer <- DESeq(dds_cancer)
results_cancer <- as.data.frame(results(dds_cancer, contrast = c("Group", "Diseased", "Healthy")))
# Replace all NA values with 1 in results_early
results_cancer[is.na(results_cancer)] <- 1
# Filter for significant genes with log2FC >= 1 or <= -1 and p-value <= 0.05
Cancer4T1_DEGs <- results_cancer[
  (results_cancer$log2FoldChange >= 1 | results_cancer$log2FoldChange <= -1) &
    results_cancer$pvalue <= 0.05,
]
write.csv(Cancer4T1_DEGs,"/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_4T1Cancer.csv")

##### EAE----
meta_EAE <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/PanDiseaseInflammation/EAE_RNASeq_Metadata.csv", sep=",", header=T) # Metadata file
meta_EAE <- as.data.frame(meta_EAE)


counts_EAE <- read.csv("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/PanDiseaseInflammation/EAE_RNASeq.csv", 
                       header = TRUE, check.names = FALSE)
counts_EAE <- na.omit(counts_EAE)

counts_EAE <- counts_EAE[!duplicated(counts_EAE[, 1]), ]
genes <- counts_EAE[, 1]
rownames(counts_EAE) <- genes
counts_EAE <- counts_EAE[, -1]

counts_EAE <- flexiDEG.function1(counts_EAE, meta_EAE, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
counts_EAE <- counts_EAE[!grepl("^Gm[0-9]", rownames(counts_EAE)), ]

# DESEQ2
dds_EAE <- DESeqDataSetFromMatrix(countData = counts_EAE, colData = meta_EAE, design = ~ Group)
dds_EAE <- DESeq(dds_EAE)
results_EAE <- as.data.frame(results(dds_EAE, contrast = c("Group", "Diseased", "Healthy")))
# Replace all NA values with 1 in results_early
results_EAE[is.na(results_EAE)] <- 1
# Filter for significant genes with log2FC >= 1 or <= -1 and p-value <= 0.05
EAE_DEGs <- results_EAE[
  (results_EAE$log2FoldChange >= 1 | results_EAE$log2FoldChange <= -1) &
    results_EAE$pvalue <= 0.05,
]
write.csv(EAE_DEGs,"/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_EAE.csv")


#### Comparison Across Multiple Disease ----

# Assuming rownames are gene names
genes_EAE <- rownames(EAE_DEGs)
genes_Cancer <- rownames(Cancer4T1_DEGs)
genes_T1D <- rownames(T1D_DEGs_Early)

# Custom colors (calmer, more refined palette)
custom_colors <- c("#3E8EDE", "#D94E5D", "#56B870")  # Blue, red, green

# Create high-quality Venn plot
venn.plot <- venn.diagram(
  x = list(
    `EAE` = genes_EAE,
    `Cancer` = genes_Cancer,
    `Type 1 Diabetes` = genes_T1D
  ),
  filename = NULL,
  fill = custom_colors,
  alpha = 0.7,
  lty = "blank", # No borders
  cex = 2.5, # Inner number font size
  fontface = "bold",
  fontfamily = "Helvetica",
  cat.fontface = "bold",
  cat.fontfamily = "Helvetica",
  cat.cex = 2.2,
  cat.col = custom_colors,
  margin = 0.08,
  main.cex = 2.5,
  main.fontface = "bold"
)

grid.newpage()
grid.draw(venn.plot)

## GSEA Analysis----


### T1D----
results_early$gene <- rownames(results_early)


gsea_NOD<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/GSEA_NODProgressorVsNonProgressor_AllTimepoints.csv",
                   row.names = 1)


library(dplyr)
library(stringr)

## 1) Prepare the NOD Early Stage table
T1D_GSEA_early <- gsea_NOD %>%
  mutate(early_Description = ifelse(
    "early_Description" %in% names(gsea_NOD) & !is.na(early_Description) & early_Description != "",
    early_Description,
    rownames(gsea_NOD)
  )) %>%
  # keep only early_ columns
  select(starts_with("early_")) %>%
  # strip the "early_" prefix so we end up with Description, NES, pvalue, etc.
  rename_with(~ sub("^early_", "", .x)) %>%
  # trim 
  mutate(Description = str_trim(Description))

# Extract results for Early Stage T1D
T1D_GSEA_early <- as.data.frame(T1D_GSEA_early)
T1D_GSEA_early$ID<-T1D_GSEA_early$Description

### Cancer----
results_cancer$gene <- rownames(results_cancer)
results_cancer <- results_cancer[order(results_cancer$stat, decreasing = TRUE), ]
lfc_vector_cancer <- setNames(results_cancer$stat, results_cancer$gene)

Cancer_GSEA <- GSEA(
  geneList = lfc_vector_cancer, 
  minGSSize = 5, 
  maxGSSize = 500, 
  pvalueCutoff = 1, 
  eps = 0,
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = mm_all_df  
)

# Extract results for Cancer
Cancer_GSEA <- as.data.frame(Cancer_GSEA)
write.csv(Cancer_GSEA, "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/Cancer_GSEA_Analysis.csv")

### EAE----
results_EAE$gene <- rownames(results_EAE)
results_EAE <- results_EAE[order(results_EAE$stat, decreasing = TRUE), ]
lfc_vector_EAE <- setNames(results_EAE$stat, results_EAE$gene)
EAE_GSEA <- GSEA(
  geneList = lfc_vector_EAE, 
  minGSSize = 5, 
  maxGSSize = 500, 
  pvalueCutoff = 1, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = mm_all_df  
)

# Extract results for EAE
EAE_GSEA <- as.data.frame(EAE_GSEA)
write.csv(EAE_GSEA, "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/EAE_GSEA_Analysis.csv")


### Combined ----

# Add disease labels
T1D_GSEA_early$Disease <- "T1D"
Cancer_GSEA$Disease <- "Cancer"
EAE_GSEA$Disease <- "EAE"

# Select relevant columns
T1D_sel <- T1D_GSEA_early %>% select(ID, NES, p.adjust, Disease)
Cancer_sel <- Cancer_GSEA %>% select(ID, NES, p.adjust, Disease)
EAE_sel <- EAE_GSEA %>% select(ID, NES, p.adjust, Disease)

# Combine into one dataframe
gsea_combined <- bind_rows(T1D_sel, Cancer_sel, EAE_sel)


filtered_pathways <- c(
  "KEGG_PATHWAYS_IN_CANCER",
  "KEGG_P53_SIGNALING_PATHWAY",
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_GLYCOLYSIS",
  "KEGG_INSULIN_SIGNALING_PATHWAY",
  "KEGG_PPAR_SIGNALING_PATHWAY",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_MTORC1_SIGNALING"
)

# Filter the combined GSEA dataframe
gsea_combined <- gsea_combined %>%
  filter(ID %in% filtered_pathways)



# Format for plotting
gsea_combined <- gsea_combined %>%
  mutate(
    neglog10padj = -log10(p.adjust ),
    ID = factor(ID, levels = rev(filtered_pathways))  
  )

# Plot
ggplot(gsea_combined, aes(x = Disease, y = ID)) +
  geom_point(aes(size = neglog10padj, color = NES)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(2,15), name = bquote(bold("-log10(" * p[adj] * ")"))) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Enrichment Across Diseases",
    x = "Disease Model",
    y = "Pathway",
    color = "NES"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14, face = "bold", angle = 30, hjust = 1),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# NOD vs AT Model----

## Adoptive Transfer - Day 10- BDC2.5 vs OVA----
#Metadata Importing
meta_AT <- read.table("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Input/metadata_AdoptiveTransfer.csv", sep=",", header=T) # Metadata file
meta_AT <- as.data.frame(meta_AT)

counts_AT <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Input/AdoptiveTransfer_Rawcounts.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_AT <- na.omit(counts_AT)

dim(counts_AT)
dim(meta_AT)
counts_AT <- flexiDEG.function1(counts_AT, meta_AT, # Run Function 1
                                convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
counts_AT <- counts_AT[!grepl("^Gm[0-9]", rownames(counts_AT)), ]

# DESeq2 analysis for Adoptive Transfer
# set references
meta_AT$Group <- relevel(factor(meta_AT$Group), ref = "OVA")
meta_AT$Time  <- relevel(factor(meta_AT$Time),  ref = "D0")

dds_AT <- DESeqDataSetFromMatrix(countData = counts_AT, colData = meta_AT,
                                 design = ~ Time + Group + Time:Group)
dds_AT <- DESeq(dds_AT)
resultsNames(dds_AT)
res_BDC_vs_OVA_d0 <- results(dds_AT, name = "Group_BDC_vs_OVA")

res_BDC_vs_OVA_d10 <- results(
  dds_AT,
  contrast = list(c("Group_BDC_vs_OVA","TimeD10.GroupBDC"))
)

AT_genes_D10  <- as.data.frame(res_BDC_vs_OVA_d10);  AT_genes_D10$gene  <- rownames(res_BDC_vs_OVA_d10)
write.csv(AT_genes_D10,
          "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_BDCvsOVA_AdoptiveTransfer_Day10.csv",
          row.names = TRUE)


# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_AT_D10  <- AT_genes_D10$stat;  names(lfc_vector_AT_D10)  <- rownames(AT_genes_D10)

# Drop NAs
lfc_vector_AT_D10  <- lfc_vector_AT_D10[!is.na(lfc_vector_AT_D10)]
lfc_vector_AT_D10  <- sort(lfc_vector_AT_D10,  decreasing = TRUE)

gsea_results_AT <- GSEA(
  geneList      = lfc_vector_AT_D10,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_AT_df <- as.data.frame(gsea_results_AT)

write.csv(gsea_results_AT_df,
          "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/GSEAResults_BDCvsOVA_AdoptiveTransfer_Day10.csv",
          row.names = FALSE)


## Double Volcano Plot- NOD-Early vs AT----

results_NOD_early<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_Early_NOD_ProgressorVsNonProgressor.csv",row.names = 1)
AT_genes_D10<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/DEG_DESEQ_BDCvsOVA_AdoptiveTransfer_Day10.csv",
                       row.names = 1)
# Replace NA p values with 1
results_NOD_early$padj[is.na(results_NOD_early$padj)] <- 1
AT_genes_D10$padj[is.na(AT_genes_D10$padj)] <- 1


# Merge the results and categorize based on thresholds
results_NOD_early$logFC_NOD <- results_NOD_early$log2FoldChange
results_NOD_early$pval_NOD <- results_NOD_early$pvalue
results_NOD_early <- results_NOD_early[, !colnames(results_NOD_early) %in% c("log2FoldChange", "pvalue")]

AT_genes_D10$logFC_AT <- AT_genes_D10$log2FoldChange
AT_genes_D10$pval_AT <- AT_genes_D10$pvalue
AT_genes_D10 <- AT_genes_D10[, !colnames(AT_genes_D10) %in% c("log2FoldChange", "pvalue")]

results_NOD_early$gene <- rownames(results_NOD_early)
AT_genes_D10$gene <- rownames(AT_genes_D10)

merged_results <- merge(results_NOD_early, AT_genes_D10, by = "gene", suffixes = c("_NOD", "_AT"))

# Thresholds
logFC_threshold <- 1
pval_threshold <- 0.05  # Adjusted p-value threshold

# Define regulation categories based on fold change and p-value cutoffs
merged_results$regulation <- "Not Significant"  # Default category
merged_results$regulation[merged_results$logFC_NOD > logFC_threshold & merged_results$pval_NOD < pval_threshold] <- "Up_NOD"
merged_results$regulation[merged_results$logFC_NOD < -logFC_threshold & merged_results$pval_NOD < pval_threshold] <- "Down_NOD"
merged_results$regulation[merged_results$logFC_AT > logFC_threshold & merged_results$pval_AT < pval_threshold] <- "Up_AT"
merged_results$regulation[merged_results$logFC_AT < -logFC_threshold & merged_results$pval_AT < pval_threshold] <- "Down_AT"
# Both Upregulated condition
merged_results$regulation[merged_results$logFC_NOD > logFC_threshold & 
                            merged_results$logFC_AT > logFC_threshold & 
                            (merged_results$pval_NOD < pval_threshold & 
                               merged_results$pval_AT < pval_threshold)] <- "Both_Up"
# Both Downregulated condition
merged_results$regulation[merged_results$logFC_NOD < -logFC_threshold & 
                            merged_results$logFC_AT < -logFC_threshold & 
                            (merged_results$pval_NOD < pval_threshold & 
                               merged_results$pval_AT < pval_threshold)] <- "Both_Down"

# Upregulated to Downregulated condition
merged_results$regulation[merged_results$logFC_NOD > logFC_threshold & 
                            merged_results$logFC_AT < -logFC_threshold & 
                            (merged_results$pval_NOD < pval_threshold & 
                               merged_results$pval_AT < pval_threshold)] <- "UpNOD_DownAT"

# Downregulated to Upregulated condition
merged_results$regulation[merged_results$logFC_NOD < -logFC_threshold & 
                            merged_results$logFC_AT > logFC_threshold & 
                            (merged_results$pval_NOD < pval_threshold & 
                               merged_results$pval_AT < pval_threshold)] <- "DownNOD_UpAT"

unique(merged_results$regulation)
write.csv(merged_results,
          "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/NODEarly_AT10_DEGs.csv",
          row.names = FALSE)



# Genes to be labelled
genes_to_label <- c( 
  "Cxcl15","Cfd","Retn",
  "Cd40lg","Cd6","Ctla4","Ms4a1","Gata1",
  "Ccr8",      
  "Fgg",
  "Ly6G",
  "Ctla4",
  "Cd163l1",    
  "Cpvl",       
  "Il23r",      
  "Lgals12",    
  "Retnlg",    
  "S100a8",     
  "S100a9",     
  "Vnn3",        
  #AT Only
  "Ctse",     
  "Cxcl13",   
  "Il31ra",   
  "Mal",      
  "Mmd2",     
  "Cd209a",   
  "Defb25",   
  "Mreg",     
  "Spib",     
  "Treml2",    
  #Both
  "Pck1",
  "Pdk4",
  #Changes
  "Fasn",
  "Lep",
  "Aldh1a7"
  
) 


ggplot(merged_results, aes(x = logFC_NOD, y = logFC_AT)) +
  # Main plot points
  geom_point(
    aes(color = regulation, size = ifelse(regulation == "Not Significant", 1, 4)),
    alpha = 0.8
  ) +
  
  # Custom color palette
  scale_color_manual(values = c(
    "Up_NOD" ="#3A6EA5",
    "Up_AT" = "#C9473A",
    "Down_NOD" = "#3A6EA5",
    "Down_AT" = "#C9473A",
    "Both_Up" = "#7F1D1D",
    "Both_Down" = "#7F1D1D",
    "UpNOD_DownAT" = "#6A3D9A",
    "DownNOD_UpAT" = "#6A3D9A",
    "Not Significant" = "gray75"
  )) +
  
  labs(
    x = expression(Log[2] ~ "FC (NOD-Progressor Vs Non-Progressor)"),
    y = expression(Log[2] ~ "FC (Adoptive Transfer-BDC2.5 vs OVA)"),
    color = "Regulation Type"
  ) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.6) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.6) +
  
  geom_label_repel(
    data = subset(merged_results, gene %in% genes_to_label),
    aes(label = gene),
    size = 6,                        
    box.padding = 0.3,
    point.padding = 0.4,               
    max.overlaps = Inf,
    force = 2,                         
    min.segment.length = 0,            
    segment.color = "black",
    segment.size = 0.8,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
    fill = "white",                    # label background
    label.size = 0.4,                  # border thickness
    label.r = unit(0.08, "lines")      # corner radius (0 for square corners)
  )+
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none",
    panel.grid = element_line(color = "grey90")
  ) +
  
  scale_size_identity() +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(-5, 5))

## GSEA-NOD Vs AT----

gsea_NOD<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/GSEA_NODProgressorVsNonProgressor_AllTimepoints.csv",
                   row.names = 1)

gsea_AT<-read.csv("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/GSEAResults_BDCvsOVA_AdoptiveTransfer_Day10.csv",
                  row.names = 1)


## Prepare the NOD Early Stage table
nod_early <- gsea_NOD %>%
  mutate(early_Description = ifelse(
    "early_Description" %in% names(gsea_NOD) & !is.na(early_Description) & early_Description != "",
    early_Description,
    rownames(gsea_NOD)
  )) %>%
  # keep only early_ columns
  dplyr::select(starts_with("early_")) %>%
  # strip the "early_" prefix so we end up with Description, NES, pvalue, etc.
  rename_with(~ sub("^early_", "", .x)) %>%
  mutate(Description = str_trim(Description))

## Prepare the AT table
at_tbl <- gsea_AT %>%
  mutate(
    Description = ifelse(
      "Description" %in% names(gsea_AT) & !is.na(Description) & Description != "",
      Description,
      rownames(gsea_AT)
    ),
    Description = str_trim(Description)
  )

## Merge on Description 
merged_NOD_AT <- inner_join(
  nod_early,
  at_tbl,
  by = "Description",
  suffix = c("_NOD", "_AT")  
)

write.csv(merged_NOD_AT,
          "/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 1/Data-Results/GSEA_NOD_AT10_Combined.csv",
          row.names = FALSE)


Immune_Metabolic_pathways <- c(
  "KEGG_PPAR_SIGNALING_PATHWAY",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_INSULIN_SIGNALING_PATHWAY",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
)


subset_df <- merged_NOD_AT %>%
  filter(Description %in% Immune_Metabolic_pathways) %>%
  dplyr::select(
    Description,
    NES_NOD, p.adjust_NOD,
    NES_AT, p.adjust_AT
  ) %>%
  as_tibble()  

plot_df <- bind_rows(
  subset_df %>%
    transmute(
      Description,
      dataset    = "NOD",
      NES        = NES_NOD,
      p.adjust   = p.adjust_NOD
    ),
  subset_df %>%
    transmute(
      Description,
      dataset    = "AT",
      NES        = NES_AT,
      p.adjust   = p.adjust_AT
    )
) %>%
  mutate(
    p.adjust = ifelse(is.na(p.adjust), 1, p.adjust),
    NES = ifelse(is.na(NES), 0, NES),
    bubble_size = -log10(p.adjust),
    ID = factor(Description, levels = Immune_Metabolic_pathways),
    dataset = factor(dataset, levels = c("NOD", "AT"))
  )



plot <- ggplot(plot_df, aes(x = dataset, y = ID)) +
  geom_point(aes(size = bubble_size, color = NES)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(5,50), name = bquote(bold("-log10(" * p[adj] * ")"))) +
  scale_y_discrete(limits = Immune_Metabolic_pathways, drop = FALSE) +  # <- enforce order
  theme_minimal() +
  labs(x = "Type 1 Diabetes Model", y = "Pathway", color = bquote(bold("NES"))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 40, face = "bold"),
    axis.text.y = element_text(size = 40),
    axis.title.x = element_text(size = 50),
    axis.title.y = element_text(size = 50),
    legend.text  = element_text(size = 30),
    legend.title = element_text(size = 40, face = "bold"),
    legend.position = "right",
    plot.title = element_text(size = 30, face = "bold")
  )

plot









