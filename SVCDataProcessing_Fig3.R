pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, dnapath, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)

library(dplyr)
# MSIGDBR Pathways ----


#Metadata Importing
meta_batch1 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Prediction Analysis/metadata_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch2 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Prediction Analysis/metadata_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch3 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Pathway Analysis/metadata_MetabolomicsCohort_PathwayAnalysis.csv", sep=",", header=T) # Metadata file
meta_batch4 <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Prediction Analysis/metadata_ScRNASeqCohort.csv", sep=",", header=T) # Metadata file
meta_test <- read.table("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Input/15039-LB_TestData/TestDataWeek6_LB_Metadata.csv", sep=",", header=T) # Metadata file

meta_batch1 <- as.data.frame(meta_batch1)
meta_batch2 <- as.data.frame(meta_batch2)
meta_batch3 <- as.data.frame(meta_batch3)
meta_batch4 <- as.data.frame(meta_batch4)
meta_test <- as.data.frame(meta_test)

# Merge metadata by columns 
meta_combined <- rbind(meta_batch1, meta_batch2,meta_batch3,meta_batch4,meta_test)
# Preview the combined metadata
head(meta_combined)

#Counts Data Importing
counts_batch1 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Prediction Analysis/gene_expected_count.annot_Jessexperimental_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch1 <- na.omit(counts_batch1)

counts_batch2 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Prediction Analysis/gene_expected_count.annot_Jessvalidation_PathwayAnalysis.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch2 <- na.omit(counts_batch2)

counts_batch3 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Prediction Analysis/gene_expected_count.annot_MetabolomicsCohort_Week6.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch3 <- na.omit(counts_batch3)

counts_batch4 <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Prediction Analysis/gene_expected_count.annot_ScRNASeqCohort.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_batch4 <- na.omit(counts_batch4)

counts_test <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Input/15039-LB_TestData/TestDataWeek6_LB.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_test <- na.omit(counts_test)


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

counts_batch4 <- counts_batch4[!duplicated(counts_batch4[, 1]), ]
genes <- counts_batch4[, 1]
rownames(counts_batch4) <- genes
counts_batch4 <- counts_batch4[, -1]


counts_test <- counts_test[!duplicated(counts_test[, 1]), ]
genes <- counts_test[, 1]
rownames(counts_test) <- genes
counts_test <- counts_test[, -1]


#Combine data
combined_counts <- merge(counts_batch1, counts_batch2, by = "row.names", all = TRUE)
combined_counts <- merge(combined_counts, counts_batch3, by.x = "Row.names", by.y = "row.names", all = TRUE)
combined_counts <- merge(combined_counts, counts_batch4, by.x = "Row.names", by.y = "row.names", all = TRUE)
combined_counts <- merge(combined_counts, counts_test, by.x = "Row.names", by.y = "row.names", all = TRUE)


# Set rownames back to genes
rownames(combined_counts) <- combined_counts$Row.names
combined_counts <- combined_counts[, -1]

# Preview the combined dataset
head(combined_counts)
# Batch Corrct Data----

early_data <- combined_counts[, meta_combined$Time == "Early"]
early_data <- na.omit(early_data)
meta_early <- meta_combined[meta_combined$Time == "Early", ]
early_data[] <- lapply(early_data, as.integer)

early_data_Batchcorrected <- flexiDEG.function1(early_data, meta_early, # Run Function 1
                                                convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                                batches = T, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0
# Remove rows where row names start with "Gm" followed by a digit
early_data_Batchcorrected <- early_data_Batchcorrected[!grepl("^Gm[0-9]", rownames(early_data_Batchcorrected)), ]

write.csv(early_data_Batchcorrected,"/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 2/Data-Input/15039-LB_TestData/Normalized_BatchCorrected_TestDataWeek6_LB.csv")




