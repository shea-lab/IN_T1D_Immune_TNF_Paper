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


# Preview the combined dataset
head(counts_AntiTNF)



# Assuming case1_f1 is your bulk expression data and meta_combined is the metadata

# Step 1: Run DESeq2 for Early and Intermediate groups
# Subset data for Early and Intermediate



counts_AntiTNF_processed <- flexiDEG.function1(counts_AntiTNF, meta_AntiTNF, # Run Function 1
                                 convert_genes = F, exclude_riken = T, exclude_pseudo = F,
                                 batches = F, quality = T, variance = F,use_pseudobulk = F) # Select filters: 0, 0,  0# Remove rows where row names start with "Gm" followed by a digit
counts_AntiTNF_processed <- counts_AntiTNF_processed[!grepl("^Gm[0-9]", rownames(counts_AntiTNF_processed)), ]


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
# Genes whose time-course differs between groups (the interaction)
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



## Venn Diagram----
deg_PreTreatment <- res_ResponderVsNonResponder_PreTreatment[
  !is.na(res_ResponderVsNonResponder_PreTreatment$pvalue) &
    res_ResponderVsNonResponder_PreTreatment$pvalue <= 0.05 &
    abs(res_ResponderVsNonResponder_PreTreatment$log2FoldChange) >= 1,
]

deg_PostTreatment <- res_ResponderVsNonResponder_PostTreatment[
  !is.na(res_ResponderVsNonResponder_PostTreatment$pvalue) &
    res_ResponderVsNonResponder_PostTreatment$pvalue <= 0.05 &
    abs(res_ResponderVsNonResponder_PostTreatment$log2FoldChange) >= 1,
]


# Gene sets
genes_PreTreatment <- rownames(deg_PreTreatment)
genes_PostTreatment <- rownames(deg_PostTreatment)

# Common + unique
genes_common      <- intersect(genes_PreTreatment, genes_PostTreatment)
genes_Pre_unique <- setdiff(genes_PreTreatment, genes_PostTreatment)
genes_Post_unique <- setdiff(genes_PostTreatment, genes_PreTreatment)

# Counts
cat("Pre-Treatment DEGs:", length(genes_PreTreatment), "\n")
cat("Post-Treatment DEGs:", length(genes_PostTreatment), "\n")
cat("Common DEGs:",     length(genes_common), "\n")
cat("Pre-unique DEGs:",length(genes_Pre_unique), "\n")
cat("Post-unique DEGs:",length(genes_Post_unique), "\n")

install.packages("ggvenn")  # run once
library(ggvenn)

venn_list <- list(
  "Pre-Treatment" = genes_PreTreatment,
  "Post-Treatment" = genes_PostTreatment
)

p_venn <- ggvenn(
  venn_list,
  fill_color = c("#B2182B", "#2166AC"),  # Nature / colorblind-safe palette
  fill_alpha = 0.4,
  stroke_size = 1.2,
  set_name_size = 6,
  text_size = 6,
  show_percentage = TRUE,
  show_elements = FALSE
) +
  theme_void() +  # removes background, axes, gridlines
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Differentially Expressed Genes\nPre vs Post Treatment")

p_venn






# Plot DEGs----


#Immune and Inflammation Related Genes
GenesPlot<-c(
             "Ly6g6c",
             "Lypd5",
             "Dsg1a",
             "Tnf",
             "Cxcl3",
             "Cxcr5",
             "Il1f5",
             "Il20ra",
             "Il22ra1",
             "Cd2",
             "Cd40lg",
             "Trat1",
             "Cd300le",
             #Down
             "Mmp7",
             "Lep",
             "Timp4",
             "Wfdc12",
             "Ces1d",
             "Retn",
            "Cd5l",
            "Klf15",
            "Epor"
)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# Thresholds
q_cut  <- 0.05
fc_cut <- 1

# Function to create color mapping
make_keyvals_fdr_fc_TNF <- function(df, q = 0.05, fc = 1,
                                          col_up = "#FF0000", col_down = "#1E90FF", col_ns = "gray70") {
  # valid stats (for coloring); everything else becomes NS
  ok   <- !is.na(df$pvalue) & !is.na(df$log2FoldChange)
  
  up   <- ok & df$pvalue < q & df$log2FoldChange >=  fc
  down <- ok & df$pvalue < q & df$log2FoldChange <= -fc
  
  # default color + label
  key   <- rep(col_ns, nrow(df))
  label <- rep("Not significant", nrow(df))
  
  # overwrite where significant
  key[up]   <- col_up
  key[down] <- col_down
  
  label[up]   <- paste0("Responder-Upregulated Post vs Pre-Treatment")
  label[down] <- paste0("Non-Responder-Upregulated Post vs Pre-Treatment")
  
  names(key) <- label        # <- legend labels; no NAs
  key
}
# keyvals_PreTreatment  <- make_keyvals_fdr_fc_TNF(PreTreatment)
# keyvals_PostTreatment <- make_keyvals_fdr_fc_TNF(PostTreatment)
keyvals_Interactions  <- make_keyvals_fdr_fc_TNF(Interactions_AntiTNF)
pick_labels <- function(df, q = 0.05, fc = 1, topN = 30) {
  idx <- which(!is.na(df$pvalue) & !is.na(df$log2FoldChange) &
                 df$pvalue < q & abs(df$log2FoldChange) >= fc)
  if (length(idx) == 0) return(character(0))
  ord <- order(df$pvalue[idx], -abs(df$log2FoldChange[idx]), na.last = NA)  # tie-break by |LFC|
  labs <- df$gene[idx][ord]
  labs <- make.unique(labs)  # avoid dup labels
  labs[seq_len(min(topN, length(labs)))]
}

# selLab_PreTreatment  <- pick_labels(PreTreatment,  q_cut, fc_cut, 30)
# selLab_PostTreatment <- pick_labels(PostTreatment, q_cut, fc_cut, 30)
#selLab_Interactions   <- pick_labels(Interactions_AntiTNF,  q_cut, fc_cut, 30)
selLab_Interactions<-GenesPlot


EnhancedVolcano(
  Interactions_AntiTNF,
  lab           = Interactions_AntiTNF$gene,
  x             = "log2FoldChange",
  y             = "pvalue",
  pCutoff       = 0.05,          # FDR threshold
  FCcutoff      = 1,
  xlab          = expression("log"[2]*"(Fold Change:Responder vs Non-Responder[Post–PreTreatment])"),
  ylab          = expression("-log"[10]*"(p)"),
  title         = "Responders vs Non-Responders(Post-Pretreatment)",
  subtitle      = paste0("p-value ≤0.05 & |LFC| ≥1 (n=", sum(Interactions_AntiTNF$pvalue<0.05 & abs(Interactions_AntiTNF$log2FoldChange)>=1, na.rm=TRUE), ")"),
  xlim          = c(-12, 12),
  ylim = c(0,6),
  boxedLabels   = TRUE,
  pointSize     = 3,
  labSize       = 6,
  colAlpha      = 0.8,
  drawConnectors= TRUE,
  max.overlaps = 50,
  colCustom     = keyvals_Interactions,
  legendPosition= "right",
  selectLab     = selLab_Interactions
)

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



# Select top DEGs(manually) related to TNF via NFKB or Macrophage/APCs

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
    aes(label = gene, color = regulation),   # <- map color here
    size = 6.5,
    box.padding = 0.3,
    point.padding = 0.4,
    max.overlaps = Inf,
    force = 2,
    min.segment.length = 0,
    segment.size = 0.6,                      # segment inherits mapped color too
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



#PLSDA Analysis----

# Prepare data
X <- t(counts_AntiTNF_processed)  # samples x genes
meta_AntiTNF <- meta_AntiTNF %>%
  mutate(SubGroup = paste(Group, Time, sep = "_"))
Y <- meta_AntiTNF$SubGroup                          # class labels
Y <- as.factor(meta_AntiTNF$SubGroup)  # Convert class labels to factor

# Now run PLS-DA
plsda_model <-  mixOmics::plsda(X, Y, ncomp = 2)
# Extract component scores
scores <- plsda_model$variates$X  
# Build data frame for ggplot
plot_df <- data.frame(LV1 = scores[,1],
                      LV2 = scores[,2],
                      Group = Y)

# Define colors & shapes
group_colors <- c(
  "Non-Responder_9" = "#B1120E",  # deep brick red
  "Non-Responder_6" = "#E74C3C",  # lighter red
  "Responder_9"     = "#1F5E4B",  # deep teal-green
  "Responder_6"     = "#6BAF92"   # lighter teal
)
group_shapes <- c("Non-Responder_6" = 16, "Non-Responder_9"=15,"Responder_6" = 16,"Responder_9"=15)


ggplot(plot_df, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.5, aes(fill = Group), show.legend = FALSE, level = 0.70) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "PLS-DA: Fast/Slow vs Non-Progressor",
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

set.seed(123)  # For reproducibility
tune_results <- tune.plsda(
  X, Y,
  ncomp = 2,
  validation = "Mfold",
  folds = 5,
  dist = "max.dist",  # or "centroids.dist"
  progressBar = TRUE
)

# Extract VIP scores for all components
vip_scores <- vip(plsda_model)
# Extract VIPs from component 1
vip_lv1 <- vip_scores[, 1]
# Export VIP scores to CSV
#write.csv(vip_scores, file = "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/Data/Sequencing/Prediction Analysis/PLSDA_Vip_scores.csv", row.names = TRUE)

# Sort descending and take top 20
top_vip_genes <- sort(vip_lv1, decreasing = TRUE)[1:20]

# Build data frame for plotting
vip_df <- data.frame(
  Gene = names(top_vip_genes),
  VIP = as.numeric(top_vip_genes)
)

# Plot
ggplot(vip_df, aes(x = reorder(Gene, VIP), y = VIP, fill = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(
    low = "#FADBD8",   # light red
    high = "#922B21"   # deep red
  ) +
  labs(
    title = "Top 100 Genes-Fast/Slow vs Non-Progressors",
    x = "Gene",
    y = "VIP Score",
    fill = "VIP"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "right"
  )
top_gene_names <- names(sort(vip_scores[,1], decreasing = TRUE)[1:20])
# Subset X to top VIP genes only
X_top <- X[, top_gene_names]
# Rerun PLS-DA
plsda_top_model <- plsda(X_top, Y, ncomp = 2)
# Extract scores
scores_top <- plsda_top_model$variates$X
# Plot
plot_df_top <- data.frame(LV1 = scores_top[,1],
                          LV2 = scores_top[,2],
                          Group = Y)

ggplot(plot_df_top, aes(x = LV1, y = LV2, color = Group, shape = Group)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Group), show.legend = FALSE, level = 0.68) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(title = "Sensitive vs Resistant (Top 20 Genes)",
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

### Heatmap Top 20----

# Expression matrix (genes x samples)
mat <- as.matrix(t(X_top))
# Metadata: group info for samples
group <- meta_AntiTNF$SubGroup
names(group) <- meta_AntiTNF$Samples

# Column annotation
ha_col <- HeatmapAnnotation(
  Group = group[colnames(mat)],
  col = list(Group = group_colors),
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  annotation_legend_param = list(title = "Group")
)

# Row-wise Z-score scaling
mat_scaled <- t(scale(t(mat)))
# Color function for heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# Heatmap with gene names
ht <- Heatmap(
  mat_scaled,
  name = "Z-score",
  top_annotation = ha_col,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 16, fontface = "italic"),  # Italicized gene names
  show_column_names = FALSE,
  heatmap_legend_param = list(title = "Expression", legend_direction = "vertical")
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()






# GSEA Analysis----
# Paths to your saved results
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
  #keyType       = "SYMBOL",       # <- tell it explicitly
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
  #keyType       = "SYMBOL",       # <- tell it explicitly
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
  #keyType       = "SYMBOL",       # <- tell it explicitly
  TERM2GENE     = mm_all_df
)
gsea_results_RvsNR_Pre_df <- as.data.frame(gsea_results_RvsNR_Pre)
write.csv(gsea_results_RvsNR_Pre_df,
          "/Users/jyotirmoyroy/Desktop/Immunometabolism T1D Paper/AntiTNFCohort/Results/GSEAResults_ResponderVsNonResponder_PreTreatment_AntiTNF.csv",
          row.names = FALSE)


library(dplyr)
library(stringr)
library(ggplot2)

# YPathways-Ineractins
Pathways_AntiTNF <- unique(c(
  "KEGG_PYRUVATE_METABOLISM",
  "KEGG_PEROXISOME",
  "KEGG_FATTY_ACID_METABOLISM",
  "KEGG_PPAR_SIGNALING_PATHWAY",
  #Down in Responders
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

# Pathways_Post_AntiTNF <- unique(c(
#   #Up
#   "KEGG_OXIDATIVE_PHOSPHORYLATION",
#   "HALLMARK_ADIPOGENESIS",
#   "HALLMARK_FATTY_ACID_METABOLISM",
#   "KEGG_LYSOSOME",
#   "KEGG_CITRATE_CYCLE_TCA_CYCLE",
#   "FAN_OVARY_CL13_MONOCYTE_MACROPHAGE",
#   "KEGG_RIBOSOME",
#   "HALLMARK_PEROXISOME",
#   "KEGG_GLYOXYLATE_AND_DICARBOXYLATE_METABOLISM",
#   "KEGG_PPAR_SIGNALING_PATHWAY",
#   "KEGG_PROPANOATE_METABOLISM",
#   "HALLMARK_BILE_ACID_METABOLISM",
#   "KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
#   "KEGG_INSULIN_SIGNALING_PATHWAY",
#   "KEGG_GLYCEROLIPID_METABOLISM",
#   "HALLMARK_PROTEIN_SECRETION",
#   "HALLMARK_XENOBIOTIC_METABOLISM",
#   "TRAVAGLINI_LUNG_MACROPHAGE_CELL",
#   "HAY_BONE_MARROW_NEUTROPHIL",
#   #Down
#   "DESCARTES_FETAL_PANCREAS_LYMPHOID_CELLS",
#   "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
#   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
#   "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
#   "HALLMARK_INFLAMMATORY_RESPONSE",
#   "HALLMARK_IL6_JAK_STAT3_SIGNALING",
#   "HE_LIM_SUN_FETAL_LUNG_C2_CXCL9_POS_MACROPHAGE_CELL",
#   "KEGG_HEMATOPOIETIC_CELL_LINEAGE",
#   "HALLMARK_IL2_STAT5_SIGNALING",
#   "AIZARANI_LIVER_C5_NK_NKT_CELLS_3",
#   "HALLMARK_KRAS_SIGNALING_UP",
#   "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
#   "TRAVAGLINI_LUNG_PLASMACYTOID_DENDRITIC_CELL",
#   "HALLMARK_APOPTOSIS",
#   "KEGG_TYPE_I_DIABETES_MELLITUS",
#   "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY"
# )
# )



# Helper to standardize clusterProfiler GSEA columns
coerce_gsea_tbl <- function(df, day_label){
  # try common column names: Description/ID/pathway/setName; p.adjust/padj
  gs  <- if ("Description" %in% names(df)) df$Description else if ("ID" %in% names(df)) df$ID else if ("pathway" %in% names(df)) df$pathway else if ("setName" %in% names(df)) df$setName else rownames(df)
  pad <- if ("p.adjust"   %in% names(df)) df$p.adjust   else if ("padj" %in% names(df)) df$padj else df$pval
  tibble(
    gs_name = as.character(gs),
    NES     = as.numeric(df$NES),
    padj    = as.numeric(pad),
    Day     = day_label
  )
}

# Build tidy tables for each day
ResponderVsNonResponder_PostVsPre_tbl  <- coerce_gsea_tbl(gsea_results_RvsNR_PostvsPre_df,  "PostVsPre")

# Keep ONLY Anti-TNF Biology Relevant
ResponderVsNonResponder_PostVsPre_tbl_subset  <- ResponderVsNonResponder_PostVsPre_tbl  %>% filter(gs_name %in% Pathways_AntiTNF)

library(dplyr)
library(ggplot2)

### Stage Wise ----

# Build tidy tables for each day
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
  # fill missing Day combos per pathway with neutral values
  tidyr::complete(gs_name, Day, fill = list(NES = 0, padj = 1)) %>%
  mutate(logp = -log10(padj))

# Order rows nicely (keep your immune_master order or order by category if you prefer)
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
  #"KEGG_CITRATE_CYCLE_TCA_CYCLE",#1 *
  "KEGG_PYRUVATE_METABOLISM",
  "KEGG_PEROXISOME",
  "KEGG_FATTY_ACID_METABOLISM",
  "KEGG_PPAR_SIGNALING_PATHWAY",
  #"KEGG_RIBOSOME",#2 *
  #Down in Responders
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_COMPLEMENT",
  #"KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",#3*
  #"KEGG_GRAFT_VERSUS_HOST_DISEASE",#4
  #"KEGG_TYPE_I_DIABETES_MELLITUS",#5
  #"HALLMARK_INTERFERON_GAMMA_RESPONSE",#6 *
  #"KEGG_CHEMOKINE_SIGNALING_PATHWAY",#7 *
  #"KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",#8 *
  #"KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",#9 *
  "KEGG_ALLOGRAFT_REJECTION",#10
  #"KEGG_JAK_STAT_SIGNALING_PATHWAY",#11 *
  #"KEGG_CELL_ADHESION_MOLECULES_CAMS",#12
  "HALLMARK_TGF_BETA_SIGNALING",
  "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  #"KEGG_TYPE_I_DIABETES_MELLITUS", #*
  "HALLMARK_IL2_STAT5_SIGNALING",
  "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
))
# *-Pre-Treatment Stage: padj >0.05- NS

# Subset to only your Anti-TNF biology pathways
AntiTNF_sets <- mm_all_sets[Pathways_AntiTNF]
#AntiTNF_sets<-mm_all_sets
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

#samples_w6 <- sample_info$Sample[sample_info$Time == "Week 6"]
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


### Plot Scores Changes Time and Group wise----
traj_df <- score_df |>
  filter(Time %in% c("Week 6", "Week 9")) |>
  mutate(
    Time = factor(Time, levels = c("Week 6","Week 9")),
    Response2 = ifelse(Response == "Responder", "Sensitive", "Resistant"),
    Response2 = factor(Response2, levels = c("Resistant","Sensitive"))
  )

paired_p <- traj_df |>
  dplyr::select(MouseID, Response2, Time, Score) |>
  tidyr::pivot_wider(names_from = Time, values_from = Score) |>
  dplyr::filter(!is.na(`Week 6`), !is.na(`Week 9`)) |>
  dplyr::group_by(Response2) |>
  dplyr::summarise(
    p = t.test(`Week 6`, `Week 9`, paired = TRUE)$p.value,
    n = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    group1  = "Week 6",
    group2  = "Week 9",
    p.label = paste0("p = ", formatC(p, format = "f", digits = 3))
  )


# choose y-positions per facet so labels don't overlap points
y_pos <- traj_df |>
  dplyr::group_by(Response2) |>
  dplyr::summarise(
    y.position = max(Score, na.rm = TRUE) +
      0.08 * diff(range(Score, na.rm = TRUE)),  # adaptive offset
    .groups = "drop"
  )

paired_p_fixed <- paired_p |>
  dplyr::select(Response2, group1, group2, p.label) |>
  dplyr::left_join(y_pos, by = "Response2")



ggplot(traj_df, aes(x = Time, y = Score, group = MouseID)) +
  geom_line(aes(color = Response2), alpha = 0.35, linewidth = 0.8) +
  geom_point(aes(color = Response2), size = 2.8) +
  facet_wrap(~Response2, nrow = 1) +
  ggpubr::stat_pvalue_manual(
    paired_p_fixed,
    label = "p.label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 6
  ) +
  scale_x_discrete(labels = c("Week 6" = "Pre", "Week 9" = "Post")) +
  scale_color_manual(values = c(
    "Resistant" = "#B2182B",
    "Sensitive" = "#3A6EA5"
  )) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Anti-TNF Response Score"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22)
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

