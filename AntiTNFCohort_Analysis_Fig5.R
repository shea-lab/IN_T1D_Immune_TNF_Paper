pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, 
               glmnet, biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, 
               apeglm, boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, 
               factoextra, edgeR, cowplot, pheatmap, coefplot, randomForest, ROCR, 
               genefilter, Hmisc, rdist, factoextra, ggforce, ggpubr, matrixStats, 
               GSEAmining, ggrepel, progress, mnormt, psych, igraph, dnapath, 
               reactome.db, GSVA, msigdbr, gglasso, MatrixGenerics, VennDiagram, 
               mikropml, glmnet, scales, stats, caret, nnet, pROC)

library(dplyr)


# Load Data ----

#Metadata Importing
meta_AntiTNF <- read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/AntiTNFCohort_metadata.csv", sep=",", header=T) # Metadata file
meta_AntiTNF <- as.data.frame(meta_AntiTNF)

#Counts Data Importing
counts_AntiTNF <- as.data.frame(read.table("/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/AntiTNFCohort_gene_counts_annot.csv", sep=",", header=T,check.names = FALSE)) # Raw counts file
counts_AntiTNF <- na.omit(counts_AntiTNF)

#Remove duplicate names
counts_AntiTNF <- counts_AntiTNF[!duplicated(counts_AntiTNF[, 1]), ]
genes <- counts_AntiTNF[, 1]
rownames(counts_AntiTNF) <- genes
counts_AntiTNF <- counts_AntiTNF[, -1]



counts_AntiTNF_processed <- flexiDEG.function1(counts_AntiTNF, meta_AntiTNF, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0# Remove rows where row names start with "Gm" followed by a digit
counts_AntiTNF_processed <- counts_AntiTNF_processed[!grepl("^Gm[0-9]", rownames(counts_AntiTNF_processed)), ]


# In Code - Responder = Sensitive ; Non-Responder = Resistant for naming schemes
# Ensure factors and correct reference levels
meta_AntiTNF$Group   <- relevel(factor(meta_AntiTNF$Group), "Non-Responder")  # reference = Responder
meta_AntiTNF$Time    <- factor(meta_AntiTNF$Time, levels = c("6", "9"))   # 6 = Pre, 9 = Post
meta_AntiTNF$MouseID <- factor(meta_AntiTNF$MouseID)
meta_AntiTNF$Batch   <- factor(meta_AntiTNF$Batch)

counts_AntiTNF_processed <- as.matrix(counts_AntiTNF_processed)
storage.mode(counts_AntiTNF_processed) <- "integer"

# Run DESEQ----

# Create DESeq2 dataset
dds_AntiTNF <- DESeqDataSetFromMatrix(
  countData = counts_AntiTNF_processed,
  colData   = meta_AntiTNF,
  design    = ~ Time+ Group+ Time:Group
)

dds_AntiTNF <- DESeq(dds_AntiTNF)
resultsNames(dds_AntiTNF)

## Interactions----
# Genes whose time-course differs between groups 
res_int <- results(dds_AntiTNF, name = "Time9.GroupResponder")
# Order by significance
res_int <- res_int[order(res_int$padj), ]
head(res_int)

## Pre-Treatment:R Vs NR----
res_ResponderVsNonResponder_PreTreatment <- results(dds_AntiTNF,
                           name = "Group_Responder_vs_Non.Responder")

## Post-Treatment:R Vs NR----
res_ResponderVsNonResponder_PostTreatment <- results(dds_AntiTNF,
                            list(c("Group_Responder_vs_Non.Responder",
                                   "Time9.GroupResponder")))

## Post Vs Pre Treatment----
res_PostVs_PreTreatment <- results(dds_AntiTNF,
                                   name = "Time_9_vs_6")
summary(res_ResponderVsNonResponder_PreTreatment)
summary(res_ResponderVsNonResponder_PostTreatment)
summary(res_int)
summary(res_PostVs_PreTreatment)


## 0) Prep results as data.frames with a 'gene' column
PreTreatment  <- as.data.frame(res_ResponderVsNonResponder_PreTreatment);  PreTreatment$gene  <- rownames(res_ResponderVsNonResponder_PreTreatment)
PostTreatment <- as.data.frame(res_ResponderVsNonResponder_PostTreatment); PostTreatment$gene <- rownames(res_ResponderVsNonResponder_PostTreatment)
Interactions_AntiTNF  <- as.data.frame(res_int);  Interactions_AntiTNF$gene  <- rownames(res_int)


write.csv(PreTreatment, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/PreTreatment_ResponderVsNonResponder_AntiTNF.csv", row.names = TRUE)
write.csv(PostTreatment, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/PostTreatment_ResponderVsNonResponder_AntiTNF.csv", row.names = TRUE)
write.csv(Interactions_AntiTNF, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/InteractionsOverTime_ResponderVsNonResponder_AntiTNF.csv", row.names = TRUE)


# Double Volcano- Pre vs Post----

RvsNR_Post_path  <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/PostTreatment_ResponderVsNonResponder_AntiTNF.csv"
RvsNR_Pre_path  <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/PreTreatment_ResponderVsNonResponder_AntiTNF.csv"

RvsNR_Post  <- read.csv(RvsNR_Post_path,  row.names = 1)
RvsNR_Pre  <- read.csv(RvsNR_Pre_path,  row.names = 1)


RvsNR_Pre$logFC_PreTreatment <- RvsNR_Pre$log2FoldChange
RvsNR_Pre$pval_PreTreatment <- RvsNR_Pre$pvalue
RvsNR_Pre <- RvsNR_Pre[, !colnames(RvsNR_Pre) %in% c("log2FoldChange", "pvalue")]

RvsNR_Post$logFC_PostTreatment <- RvsNR_Post$log2FoldChange
RvsNR_Post$pval_PostTreatment <- RvsNR_Post$pvalue
RvsNR_Post <- RvsNR_Post[, !colnames(RvsNR_Post) %in% c("log2FoldChange", "pvalue")]

RvsNR_Pre$gene <- rownames(RvsNR_Pre)
RvsNR_Post$gene <- rownames(RvsNR_Post)

# Merge the results
merged_results <- merge(RvsNR_Pre, RvsNR_Post, by = "gene", suffixes = c("_PreTreatment", "_PostTreatment"))

# Thresholds
logFC_threshold <- 1
pval_threshold <- 0.05  # p-value threshold

# Define regulation categories based on fold change and p-value cutoffs
merged_results$regulation <- "Not Significant"  # Default category
merged_results$regulation[merged_results$logFC_PreTreatment > logFC_threshold & merged_results$pval_PreTreatment < pval_threshold] <- "Up_PreTreatment"
merged_results$regulation[merged_results$logFC_PreTreatment < -logFC_threshold & merged_results$pval_PreTreatment < pval_threshold] <- "Down_PreTreatment"
merged_results$regulation[merged_results$logFC_PostTreatment > logFC_threshold & merged_results$pval_PostTreatment < pval_threshold] <- "Up_PostTreatment"
merged_results$regulation[merged_results$logFC_PostTreatment < -logFC_threshold & merged_results$pval_PostTreatment < pval_threshold] <- "Down_PostTreatment"
# Both Upregulated condition
merged_results$regulation[merged_results$logFC_PreTreatment > logFC_threshold & 
                            merged_results$logFC_PostTreatment > logFC_threshold & 
                            (merged_results$pval_PreTreatment < pval_threshold & 
                               merged_results$pval_PostTreatment < pval_threshold)] <- "Both_Up"
# Both Downregulated condition
merged_results$regulation[merged_results$logFC_PreTreatment < -logFC_threshold & 
                            merged_results$logFC_PostTreatment < -logFC_threshold & 
                            (merged_results$pval_PreTreatment < pval_threshold & 
                               merged_results$pval_PostTreatment < pval_threshold)] <- "Both_Down"

# Upregulated to Downregulated condition
merged_results$regulation[merged_results$logFC_PreTreatment > logFC_threshold & 
                            merged_results$logFC_PostTreatment < -logFC_threshold & 
                            (merged_results$pval_PreTreatment < pval_threshold & 
                               merged_results$pval_PostTreatment < pval_threshold)] <- "UpPre_DownPost"

# Downregulated to Upregulated condition
merged_results$regulation[merged_results$logFC_PreTreatment < -logFC_threshold & 
                            merged_results$logFC_PostTreatment > logFC_threshold & 
                            (merged_results$pval_PreTreatment < pval_threshold & 
                               merged_results$pval_PostTreatment < pval_threshold)] <- "DownPre_UpPost"

unique(merged_results$regulation)



# Select top DEGs related to TNF via NFKB 

selected_genes<-c( #Pre-Treatment 
  "Retn", "Retnla", "Lcn2", "Cd5l", "Chit1", "Olr1", "Lgals12", "Cfd", "Apol6", "Vnn1", "Vnn3", "Serpina3c", "Mmp10", "Cemip", "Pglyrp4", "Il1f5", "Il17b", "Il22ra1", "Adora1", "P2ry4", "Agt",
  #Post-Treatment 
  "Ido1", "Ido2", "Sectm1a", "Sectm1b", #"Pglyrp4", 
  "Serpina3g", #"Cemip", 
  "Grem1", "Il17re", "Smpd3", "Lipg", "Duox2", "Trim30c", "Slamf1", "Cd209c", "Cd300e", "Lag3", "Gpr18", "Gbp8", "Gbp2b", "Ifit1bl2", "Inpp5j",
  #Both
  "Chodl", "Mrgprg" )


top_genes <- merged_results %>%
  filter(
    gene %in% selected_genes
  ) %>%
  mutate(
    max_abs_logFC = pmax(abs(logFC_PreTreatment), abs(logFC_PostTreatment))
  ) %>%
  arrange(desc(max_abs_logFC))

ggplot(merged_results, aes(x = logFC_PreTreatment, y = logFC_PostTreatment)) +
  geom_point(
    aes(color = regulation, size = ifelse(regulation == "Not Significant", 1, 3.5)),
    alpha = 0.9
  ) +
  scale_color_manual(values = c(
    "Up_PreTreatment"   = "#3A6EA5",
    "Up_PostTreatment"  = "#C9473A",
    "Down_PreTreatment" = "#3A6EA5",
    "Down_PostTreatment"= "#C9473A",
    "Both_Up"           = "#7F1D1D",
    "Both_Down"         = "#7F1D1D",
    "UpPre_DownPost"    = "#6A3D9A",
    "DownPre_UpPost"    = "#6A3D9A",
    "Not Significant"   = "gray75"
  )) +
  labs(
    x = expression(Log[2] ~ "Fold Change (Pre-Treatment)"),
    y = expression(Log[2] ~ "Fold Change (Post-Treatment)"),
    color = "Regulation Type"
  ) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = c(-1.0, 1.0), linetype = "dashed", color = "black", linewidth = 0.6) +
  ggrepel::geom_label_repel(
    data = top_genes,
    aes(label = gene, color = regulation),   
    size = 6.5,
    box.padding = 0.3,
    point.padding = 0.4,
    max.overlaps = Inf,
    force = 2,
    min.segment.length = 0,
    segment.size = 0.6,                     
    fill = "white",
    label.size = 0.25,
    label.r = unit(0.15, "lines"),
    arrow = grid::arrow(length = unit(0.02, "npc"), type = "closed")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    legend.position = "right",
    panel.grid = element_line(color = "grey90")
  ) +
  scale_size_identity() +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(-5, 5))





# GSEA Analysis----

RvsNR_PostvsPre_path  <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/InteractionsOverTime_ResponderVsNonResponder_AntiTNF.csv"
RvsNR_Post_path  <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/PostTreatment_ResponderVsNonResponder_AntiTNF.csv"
RvsNR_Pre_path  <- "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/PreTreatment_ResponderVsNonResponder_AntiTNF.csv"

# Import
RvsNR_PostvsPre  <- read.csv(RvsNR_PostvsPre_path,  row.names = 1)
RvsNR_Post  <- read.csv(RvsNR_Post_path,  row.names = 1)
RvsNR_Pre  <- read.csv(RvsNR_Pre_path,  row.names = 1)


# Spot-check a known rejection gene (e.g., Ifng) if present:
RvsNR_PostvsPre["Ifng", c("log2FoldChange","stat","padj")]


# Build ranked gene lists using DESeq2 Wald stat
lfc_vector_RvsNR_PostvsPre  <- RvsNR_PostvsPre$stat;  names(lfc_vector_RvsNR_PostvsPre)  <- rownames(RvsNR_PostvsPre)
lfc_vector_RvsNR_Post  <- RvsNR_Post$stat;  names(lfc_vector_RvsNR_Post)  <- rownames(RvsNR_Post)
lfc_vector_RvsNR_Pre  <- RvsNR_Pre$stat;  names(lfc_vector_RvsNR_Pre)  <- rownames(RvsNR_Pre)


# Drop NAs
lfc_vector_RvsNR_PostvsPre  <- lfc_vector_RvsNR_PostvsPre[!is.na(lfc_vector_RvsNR_PostvsPre)]
lfc_vector_RvsNR_Post <- lfc_vector_RvsNR_Post[!is.na(lfc_vector_RvsNR_Post)]
lfc_vector_RvsNR_Pre <- lfc_vector_RvsNR_Pre[!is.na(lfc_vector_RvsNR_Pre)]


# Sort decreasing (required by clusterProfiler::GSEA)
lfc_vector_RvsNR_PostvsPre  <- sort(lfc_vector_RvsNR_PostvsPre,  decreasing = TRUE)
lfc_vector_RvsNR_Post <- sort(lfc_vector_RvsNR_Post,  decreasing = TRUE)
lfc_vector_RvsNR_Pre <- sort(lfc_vector_RvsNR_Pre,  decreasing = TRUE)



library(msigdbr)
library(dplyr)

# --- Collect each set and convert into 2-column (gs_name, gene_symbol) ---

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

# --- Final combined TERM2GENE data frame ---
mm_all_df <- rbind(mm_h_df, mm_kegg_df)

library(clusterProfiler)
library(msigdbr)

gsea_results_RvsNR_PostvsPre <- GSEA(
  geneList      = lfc_vector_RvsNR_PostvsPre,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_RvsNR_PostvsPre_df <- as.data.frame(gsea_results_RvsNR_PostvsPre)
write.csv(gsea_results_RvsNR_PostvsPre_df,
          "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/GSEAResults_ResponderVsNonResponder_PostVsPreTreatment_AntiTNF.csv",
          row.names = FALSE)


gsea_results_RvsNR_Post <- GSEA(
  geneList      = lfc_vector_RvsNR_Post,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_RvsNR_Post_df <- as.data.frame(gsea_results_RvsNR_Post)
write.csv(gsea_results_RvsNR_Post_df,
          "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/GSEAResults_ResponderVsNonResponder_PostTreatment_AntiTNF.csv",
          row.names = FALSE)


gsea_results_RvsNR_Pre <- GSEA(
  geneList      = lfc_vector_RvsNR_Pre,
  minGSSize     = 5,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  eps           = 0,
  seed          = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE     = mm_all_df
)
gsea_results_RvsNR_Pre_df <- as.data.frame(gsea_results_RvsNR_Pre)
write.csv(gsea_results_RvsNR_Pre_df,
          "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/GSEAResults_ResponderVsNonResponder_PreTreatment_AntiTNF.csv",
          row.names = FALSE)


library(dplyr)
library(stringr)
library(ggplot2)

# Selected Pahways to plot
Pathways_AntiTNF <- unique(c(
  "KEGG_PYRUVATE_METABOLISM",
  "KEGG_PEROXISOME",
  "KEGG_FATTY_ACID_METABOLISM",
  "KEGG_PPAR_SIGNALING_PATHWAY",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_TGF_BETA_SIGNALING",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_TYPE_I_DIABETES_MELLITUS",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  ))


coerce_gsea_tbl <- function(df, day_label){
  gs  <- if ("Description" %in% names(df)) df$Description else if ("ID" %in% names(df)) df$ID else if ("pathway" %in% names(df)) df$pathway else if ("setName" %in% names(df)) df$setName else rownames(df)
  pad <- if ("p.adjust"   %in% names(df)) df$p.adjust   else if ("padj" %in% names(df)) df$padj else df$pval
  tibble(
    gs_name = as.character(gs),
    NES     = as.numeric(df$NES),
    padj    = as.numeric(pad),
    Day     = day_label
  )
}




### Stage Wise Plot----

# Build tidy tables for each tim poiny
PreTreatment_tbl  <- coerce_gsea_tbl(gsea_results_RvsNR_Pre_df,  "Pre-Treatment")
PostTreatment_tbl <- coerce_gsea_tbl(gsea_results_RvsNR_Post_df, "Post-Treatment")

# Keep ONLY Pathways_AntiTNF pathways
PreTreatment_sel  <- PreTreatment_tbl  %>% filter(gs_name %in% Pathways_AntiTNF)
PostTreatment_sel <- PostTreatment_tbl %>% filter(gs_name %in% Pathways_AntiTNF)

# Combine and ensure both days present for every pathway
plot_df <- bind_rows(PreTreatment_sel, PostTreatment_sel) %>%
  mutate(
    gs_name = as.character(gs_name),
    Day = factor(Day, levels = c("Pre-Treatment","Post-Treatment"))
  ) %>%
  tidyr::complete(gs_name, Day, fill = list(NES = 0, padj = 1)) %>%
  mutate(logp = -log10(padj))


plot_df$gs_name <- factor(plot_df$gs_name, levels = rev(Pathways_AntiTNF))

# Plot: rows = pathways, columns = Day; size = −log10(padj); color = NES
p <- ggplot(plot_df, aes(x = Day, y = gs_name)) +
  geom_point(aes(size = logp, color = NES)) +
  scale_size_continuous(name = "−log10(padj)", range = c(2, 10)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  labs(x = "", y = "", title = "Anti-TNF Treated: Sensitive Vs Resistant") +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  )
print(p)



# GSVA-TNF Scoring----
library(dplyr)
library(ggplot2)
library(ggpubr)

mm_all_sets <- c(mm_h_sets, mm_kegg_sets)
Pathways_AntiTNF <- unique(c(
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_PYRUVATE_METABOLISM",
  "KEGG_PEROXISOME",
  "KEGG_FATTY_ACID_METABOLISM",
  "KEGG_PPAR_SIGNALING_PATHWAY",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "KEGG_ALLOGRAFT_REJECTION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
))


# Subset to only TNF related biology pathways
AntiTNF_sets <- mm_all_sets[Pathways_AntiTNF]

AntiTNF_sets <- AntiTNF_sets[!vapply(AntiTNF_sets, is.null, logical(1))]
vst_obj  <- vst(dds_AntiTNF, blind = FALSE)
expr_mat <- assay(vst_obj)

gsvapar <- gsvaParam(
  expr     = expr_mat,
  geneSets = AntiTNF_sets,
  maxDiff  = TRUE
)

gsva_AntiTNF <- gsva(gsvapar)

sample_info <- colData(dds_AntiTNF) |> as.data.frame()

sample_info <- sample_info |>
  transmute(
    Sample   = Samples,
    MouseID  = MouseID,
    Time     = ifelse(Time == 6, "Week 6", "Week 9"),
    Response = factor(Group, levels = c("Non-Responder", "Responder"))
  )

rownames(sample_info) <- sample_info$Sample

samples_w6 <- sample_info$Sample[sample_info$Time %in% c("Week 6")]

X6 <- t(gsva_AntiTNF[, samples_w6])
y6 <- sample_info[samples_w6, "Response"]
y6_bin <- as.numeric(y6 == "Responder")

set.seed(123)
foldid <- seq_along(y6_bin)  # LOOCV


cvfit <- cv.glmnet(
  x = X6,
  y = y6_bin,
  family = "binomial",
  alpha = 0.7,
  foldid = foldid,
  standardize = TRUE
)

coef_final <- coef(cvfit, s = "lambda.min")

selected_pathways <- setdiff(
  rownames(coef_final)[coef_final[,1] != 0],
  "(Intercept)"
)
selected_pathways
beta      <- coef_final[selected_pathways, 1]
intercept <- coef_final["(Intercept)", 1]

score_samples <- function(samples) {
  X <- t(gsva_AntiTNF[selected_pathways, samples, drop = FALSE])
  as.numeric(intercept + X %*% beta)
}

score_df <- sample_info |>
  mutate(
    Score = score_samples(Sample)
  )

score_df <- score_df |>
  group_by(Time) |>
  mutate(Score_Z = as.numeric(scale(Score))) |>
  ungroup()

### Plot Baseline Scores----

w6_df <- score_df |>
  filter(Time == "Week 6") |>
  mutate(Response2 = ifelse(Response == "Responder", "Sensitive", "Resistant"),
         Response2 = factor(Response2, levels = c("Resistant", "Sensitive")))

p_w6 <- t.test(Score ~ Response2, data = w6_df)$p.value

p_df <- data.frame(
  group1 = "Resistant",
  group2 = "Sensitive",
  y.position = max(w6_df$Score, na.rm = TRUE) +
    0.08 * diff(range(w6_df$Score, na.rm = TRUE)),
  p.label = paste0("p = ", formatC(p_w6, format = "f", digits = 3))
)

ggplot(w6_df, aes(x = Response2, y = Score, fill = Response2)) +
  geom_boxplot(width = 0.55, alpha = 0.75,
               outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.12, size = 2.8, alpha = 0.9) +
  ggpubr::stat_pvalue_manual(
    p_df,
    label = "p.label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    size = 7
  ) +
  scale_fill_manual(values = c(
    "Resistant" = "#B2182B",
    "Sensitive" = "#3A6EA5"
  )) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Anti-TNF Response Score",
    title = "Pre-treatment"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    plot.title = element_text(size = 22, face = "bold")
  )




### Plot Delta Scores----
delta_df <- traj_df |>
  dplyr::select(MouseID, Response2, Time, Score) |>
  tidyr::pivot_wider(names_from = Time, values_from = Score) |>
  mutate(Delta_ARPS = `Week 9` - `Week 6`)

ggplot(delta_df, aes(x = Response2, y = Delta_ARPS, fill = Response2)) +
  geom_boxplot(width = 0.55, alpha = 0.75, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.12, size = 2.8, alpha = 0.9) +
  ggpubr::stat_compare_means(
    comparisons = list(c("Resistant","Sensitive")),
    method = "t.test",
    label = "p.format",   # << numeric p-value
    digits = 3,           # << exactly 3 decimals
    tip.length = 0.02,
    size = 7
  ) +
  scale_fill_manual(values = c("Resistant"="#B2182B","Sensitive"="#3A6EA5")) +
  theme_classic(base_size = 16) +
  labs(x = NULL, y = expression(Delta*"Anti-TNF Response Score"),
       title = "Response-Time interaction") +
  theme(legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 18, face = "bold"))

write.csv(score_df,"/Users/jyotirmoyroy/Desktop/T1S_ImmunometabolismPaper/Figure 5/Data-Results/ATRSScore.csv")

