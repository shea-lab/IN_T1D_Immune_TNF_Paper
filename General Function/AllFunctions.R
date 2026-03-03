# Gene conversion functions ---- 
# These functions are called out in the first central function of the package
convMensID.Msym <- function(counts, mouse_mart = NULL) { # Convert mouse IDs to symbols
  genes <- rownames(counts) # counts df must have genes as rownames, not the first column
  require("biomaRt") # Necessary package
  if (is.null(mouse_mart)) { # Default mart
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE,
                     host="https://dec2021.archive.ensembl.org/")
  } else {
    mouse <- mouse_mart # Mart specified by user
  }
  new_genes <- getLDS(attributesL = c("mgi_symbol"), filters = "ensembl_gene_id", 
                      values = genes, mart = mouse, attributes = c("ensembl_gene_id"), 
                      martL = mouse, uniqueRows=TRUE)
  colnames(new_genes) <- c("old", "new")
  new_genes <- new_genes[!duplicated(new_genes[1]),]
  new_genes <- new_genes[!duplicated(new_genes[2]),]
  rownames(new_genes) <- new_genes$old
  counts_conv <- merge(new_genes, counts, by = "row.names", all.x=T)
  rownames(counts_conv) <- counts_conv$new
  counts_conv <- counts_conv[,-1:-3]
  return(counts_conv)
}
convHensID.Hsym <- function(counts, human_mart = NULL) { # Convert human IDs to symbols
  genes <- rownames(counts) # counts df must have genes as rownames, not the first column
  require("biomaRt") # Necessary package
  if (is.null(human_mart)) { # Default mart
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE,
                     host="https://dec2021.archive.ensembl.org/")
  } else {
    human <- human_mart # Mart specified by user
  }
  new_genes <- getLDS(attributesL = c("hgnc_symbol"), filters = "ensembl_gene_id", 
                      values = genes, mart = human, attributes = c("ensembl_gene_id"), 
                      martL = human, uniqueRows=TRUE)
  colnames(new_genes) <- c("old", "new")
  new_genes <- new_genes[!duplicated(new_genes[1]),]
  new_genes <- new_genes[!duplicated(new_genes[2]),]
  rownames(new_genes) <- new_genes$old
  counts_conv <- merge(new_genes, counts, by = "row.names", all.x=T)
  rownames(counts_conv) <- counts_conv$new
  counts_conv <- counts_conv[,-1:-3]
  return(counts_conv)
}

# Central functions ----

# In the first function, counts are filtered as specified.
flexiDEG.function1 <- function(counts, # Counts df, genes=rows & samples=cols, 
                               # gene names can be row names or the first col
                               metadata, # Metadata with first cell in the col reading "Samples"
                               # also needs a col with first cell reading "Group"
                               # and a col with first cell reading "Batch" if going to batch correct
                               convert_genes = FALSE, # Convert gene IDs to symbols 
                               exclude_riken = TRUE, # Remove riken genes prior to analyses
                               exclude_pseudo = TRUE, # Remove psuedo genes prior to analyses 
                               batches = FALSE, # Analyze multiple data sets
                               quality = FALSE, # Report on data quality
                               variance = TRUE, # Filter out genes that change little
                               use_pseudobulk = FALSE # New argument to specify pseudobulk usage
) {
  # Check inputs ---- 
  sample_names <- metadata$Samples # Extract sample names 
  sample_num <- length(sample_names) # Extract number of samples 
  groups <- metadata$Group # Extract sample groups 
  message <- "Please indicate if the genes are rownames (1) or the first column (2) of the counts df: "
  user_input <- readline(prompt = message) # Check for consistency between metadata & counts
  if (user_input == 1) {
    if (ncol(counts) != nrow(metadata)) {stop("Number of samples in counts & metadata are inconsistent")}
    if (!all(sample_names %in% colnames(counts))) {stop("Sample names in counts & metadata are inconsistent")} 
    genes <- make.names(rownames(counts), unique = T) # Cant start rows w/ numbers, add a letter to genes starting w/ numbers
  } else if (user_input == 2) {
    if (ncol(counts) - 1 != nrow(metadata)){stop("Number of samples in counts & metadata are inconsistent")}
    if (!all(sample_names %in% colnames(counts[, -c(1)]))) {stop("Sample names in counts & metadata are inconsistent")} 
    genes <- make.names(counts[, 1], unique = T) # Cant start rows w/ numbers, add a letter to genes starting w/ numbers
    counts <- counts[, -1] # Remove first column of genes to add it back as row names instead 
    colnames(counts) <- c(sample_names) # Add sample names as column names
    rownames(counts) <- genes # Add gene names as row names
  } else {stop("Genes must be listed as rownames or first column in the counts")
  }
  if (sum(duplicated(genes)) != 0) { # Check for duplicate genes
    genes2 <- genes[!duplicated(genes)] # Genes which are not duplicated
    counts <- counts[genes2, ] # Remove duplicated genes
    removed <- length(genes) - length(genes2) # Number of removed genes
    print(paste(removed, " duplicate genes removed")) # Display number 
    genes <- rownames(counts) # Redefine genes without duplicates
  } 
  
  # Convert genes ---- 
  if (convert_genes == TRUE) { # Only converts ensembl IDs to symbols for mouse or human
    if (substring(genes[1], 1, 3) != "ENS") { # Check for genes as ensembl IDs
      stop("Only genes as ensembl IDs can be converted by this function.")
    } else {
      if (substring(genes[1], 4, 7) == "MUSG") { # If they are mouse
        counts <- convMensID.Msym(counts) # Function to convert
        dt <- "mouse"
      } else if (substring(genes[1], 4, 4) == "G") { # If they are human
        counts <- convHensID.Hsym(counts) # Function to convert
        dt <- "human"
      } else { 
        stop("Genes are not mouse or human; convert manually and then input into function.")
      }}
    genes <- rownames(counts) # Redefine genes once converted
  } else { # Check if genes are human or mouse, even though we won't be converting
    letters_before_numbers <- gsub("[^A-Za-z]+.*$", "", genes[1:10]) # Find letters before the first number
    all_caps <- all(grepl("^[A-Z]+$", letters_before_numbers)) # Check if all characters are uppercase
    if (all_caps == TRUE) {
      dt <- "human"
    } else {
      dt <- "mouse"
    }}
  
  # Exclude genes ---- 
  if (exclude_riken == TRUE | exclude_pseudo == TRUE) { 
    if (substring(genes[1], 1, 3) == "ENS") { # Check for genes as ensembl IDs
      stop("Riken genes cannot be excluded when given as ensembl IDs")
    }}
  if (exclude_riken == TRUE) { 
    if (dt == "mouse") {
      genes2 <- genes[grep("Rik", genes, invert = T)] # Find Riken genes
      removed <- length(genes) - length(genes2) # Number of removed genes
      print(paste(removed, "Riken genes removed")) # Display number 
      counts <- counts[genes2,] # Remove excluded genes
      genes <- rownames(counts)  # Redefine genes without Riken genes
    } else if (dt == "human") {
      genes2 <- genes[grep("Rik", genes, invert = T)] # Find Riken genes
      removed <- length(genes) - length(genes2) # Number of removed genes
      print(paste(removed, "RIKEN genes removed")) # Display number 
      counts <- counts[genes2,] # Remove excluded genes
      genes <- rownames(counts)  # Redefine genes without Riken genes
    }}
  if (exclude_pseudo == TRUE) { 
    if (dt == "mouse") {
      genes2 <- genes[grep("Ps", genes, invert = T)] # Find Pseudo genes
      genes2 <- genes2[grep("Gm", genes2, invert = T)] # Find Pseudo genes
      removed <- length(genes) - length(genes2)  # Number of removed genes
      print(paste(removed, "Pseudo genes (starting with Ps or Gm) removed")) 
      counts <- counts[genes2,] # Remove excluded genes
      genes <- rownames(counts)  # Redefine genes without Pseudo genes
    } else if (dt == "human") {
      genes2 <- genes[grep("PS", genes, invert = T)] # Find Pseudo genes
      genes2 <- genes2[grep("GM", genes2, invert = T)] # Find Pseudo genes
      removed <- length(genes) - length(genes2)  # Number of removed genes
      print(paste(removed, "Pseudo genes (starting with PS or GM) removed")) 
      counts <- counts[genes2,] # Remove excluded genes
      genes <- rownames(counts)  # Redefine genes without Pseudo genes
    }}
  
  if (use_pseudobulk== FALSE){
    # Remove unexpressed genes ---- 
    counts[is.na(counts)] <- 0 # Set any missing values equal to zero
  }
  else{
    # Check fraction of missing values in counts ----
    missing_fraction <- mean(is.na(counts)) # Calculate fraction of missing values
    print(paste("Fraction of missing values in counts data:", round(missing_fraction, 4)))
    
    # Function to filter genes based on a missing value threshold
    filter_genes <- function(df, max_missing_threshold = 0.2) {
      keep_genes <- apply(df, 1, function(x) mean(is.na(x)) <= max_missing_threshold)
      df_filtered <- df[keep_genes, ]
      return(df_filtered)
    }
    message <- "Please maximum fraction of missing values allowed: "
    max_missing_threshold <- as.numeric(readline(prompt = message))
    counts  <- filter_genes(counts , max_missing_threshold = 0.2)
    # Impute missing values with median across samples
    median_impute <- function(df) {
      df[] <- lapply(df, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
      return(df)
    }
    counts <- median_impute(counts)
  }
  counts2 <- counts[rowSums(counts) > 0,] # Remove genes unexpressed across all samples
  removed <- nrow(counts) - nrow(counts2) # Number of removed genes
  print(paste(removed, "unexpressed genes removed")) # Display number
  counts <- counts2 # Update counts df
  counts2 <- counts[rowSums(counts == 0) <= sample_num*(0.85),] # Unexpressed in 85%+ of samples
  removed <- nrow(counts) - nrow(counts2) # Number of removed genes
  print(paste(removed, "mostly zero genes removed")) # Display number
  counts <- counts2 # Update counts df
  
  # Batch correction ---- 
  if (batches == TRUE) { # Batch correction requires three things:
    # 1) counts df w/ RAW data from all batches combined, rows=genes cols=samples
    # 2) vector defining sample batches
    # 3) vector defining sample groups
    Batch <- metadata$Batch 
    coldata <- data.frame(cbind(groups, Batch)) # Combine groups and batch 
    row.names(coldata) <- sample_names # DESeq2 uses raw counts; rows:genes & cols:samples
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Batch)
    paste(nrow(dds), " genes input for batch correction", sep="")
    suppressMessages(dds <- DESeq(dds)) # Analyze previously specified variables
    vsd <- vst(dds, blind = FALSE) # Estimates dispersion trend & stabilizes transformation
    mat <- assay(vsd)
    mat <- limma::removeBatchEffect(mat, vsd$Batch)
    assay(vsd) <- mat
    counts <- as.data.frame(assay(vsd)) # Summarized assay value, Update counts df
  }
  
  #ps_counts <- log(counts + 1, base = 2) # Pseudo-normalize
  ps_counts<- counts
  # All filters are per sample per gene
  message <- "Please enter expression cutoff for low filter (required): "
  filter_low <- as.numeric(readline(prompt = message))
  if (is.na(filter_low)) {stop("Invalid input. Please enter a numeric value.")} 
  keep1 <- rowSums(ps_counts) > (filter_low*sample_num) # Low count filter
  ps_filt <- ps_counts[keep1,] # Apply filter to pseudonorm counts
  removed <- nrow(ps_counts) - nrow(ps_filt) # Number of filtered genes
  print(paste(removed, "genes removed due to expression below the low filter"))
  message <- "Please enter expression cutoff for group filter (0 for no group filter): "
  filter_group <- as.numeric(readline(prompt = message))
  if (is.na(filter_group)) {stop("Invalid input. Please enter a numeric value.")
  } else if (filter_group != 0) { 
    keep2 <- rowSums(ps_filt) >= (filter_group*sample_num) # Group filter, optional
    ps_filt2 <- ps_filt[keep2,] # Apply filter to pseudonorm counts
    removed <- nrow(ps_filt) - nrow(ps_filt2) # Number of filtered genes
    print(paste(removed, "genes removed due to expression below the low group filter"))
    ps_filt <- ps_filt2 # Update counts df
  }
  message <- "Please enter expression cutoff for high filter (0 for no high filter): " 
  filter_high <- as.numeric(readline(prompt = message))
  if (is.na(filter_high)) {stop("Invalid input. Please enter a numeric value.")
  } else if (filter_high != 0) { 
    keep3 <- rowSums(ps_filt) < (filter_high*sample_num) # High count filter, optional
    ps_filt2 <- ps_filt[keep3,] # Apply filter to pseudonorm counts
    removed <- nrow(ps_filt) - nrow(ps_filt2) # Number of filtered genes
    print(paste(removed, "genes removed due to expression above the high filter"))
    ps_filt <- ps_filt2 # Update counts df
  }
  genes_filt <- rownames(ps_filt) # Update genes vector for downstream 
  counts_filt <- counts[genes_filt,] # Update counts df for downstream 
  
  # Sample quality ---- 
  if (quality == TRUE) { 
    suppressMessages(ps_LF <- reshape2::melt(ps_counts, variable.name = "Samples", value.name = "Count")) # Long form
    ps_LF <- merge(ps_LF, metadata, by = "Samples") # Merge w/ metadata to get the group for each sample
    ps_LF <- ps_LF[!is.na(ps_LF$Count), ] # Remove any missing values
    suppressMessages(psf_LF <- reshape2::melt(ps_filt, variable.name = "Samples", value.name = "Count")) # Long form
    psf_LF <- merge(psf_LF, metadata, by = "Samples") # Merge w/ metadata to get the group for each sample
    psf_LF <- psf_LF[!is.na(psf_LF$Count), ] # Remove any missing values
    if (sample_num <= 20) { # If 20 samples or less, print to 1 page
      par(mfrow=c(4,ceiling(sample_num/4))) # All histograms on a page w/ 4 rows
    } else {par(mfrow=c(3,5))} # Histograms on multiple pages w/ 3 rows & 5 cols
    par(mar = c(4,1,2,1), oma = c(0.5,1,0,0)) # Adjust margins to better use space
    sapply(1:sample_num, function(x) hist(counts[,x], ylim=c(0,1000), breaks=100,
                                          xlab=sample_names[x], ylab = "", main = ""))
    mtext("Counts, 100 bins", side = 3, line = -2, outer = T) # Add title
    sapply(1:sample_num, function(x) hist(ps_counts[,x], ylim=c(0,1000), breaks=100,
                                          xlab=sample_names[x], ylab = "", main = ""))
    mtext("Log Normalized Counts, 100 bins", side = 3, line = -2, outer = T) # Add title
    sapply(1:sample_num, function(x) hist(ps_filt[,x], ylim=c(0,1000), breaks=100,
                                          xlab=sample_names[x], ylab = "", main = ""))
    mtext("Log Normalized Counts Filt, 100 bins", side = 3, line = -2, outer = T) # Add title
    # Box & Density Plots, how samples compare in mean, variance & any outliers
    p1 <- ggplot(ps_LF, aes(x = Samples, y = Count)) + 
      geom_boxplot() + xlab("") +  ylab("Log Normalized Counts") +
      scale_x_discrete(breaks = sample_names, labels = sample_names) + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
      ggtitle("Counts Distribution")
    grid.draw(p1) # Plot distribution
    p2 <- ggplot(ps_LF, aes(x = Count, fill = Group, group = Group)) + 
      geom_density(alpha = 0.4, color = "black") + 
      theme(legend.position = "top") + 
      xlab("Log Normalized Counts") + 
      ggtitle("Combined Density Distribution of Counts")
    grid.draw(p2) # Plot distribution
    # Box & Density Plots, how samples compare in mean, variance & any outliers
    p1 <- ggplot(psf_LF, aes(x = Samples, y = Count)) + 
      geom_boxplot() + xlab("") +  ylab("Log Normalized Counts") +
      scale_x_discrete(breaks = sample_names, labels = sample_names) + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
      ggtitle("Counts Distribution")
    grid.draw(p1) # Plot distribution
    p2 <- ggplot(psf_LF, aes(x = Count, fill = Group, group = Group)) + 
      geom_density(alpha = 0.4, color = "black") + 
      theme(legend.position = "top") + 
      xlab("Log Normalized Counts") + 
      ggtitle("Combined Density Distribution of Counts")
    grid.draw(p2) # Plot distribution
  }
  
  # Variance filtering ---- 
  if (variance == TRUE) { 
    # Filter by variance across all samples to exclude genes that don't change much
    mn_cut <- 0.2 # Mean cutoff
    var_cut <- 0.01 # Variance cutoff
    df <- cbind.data.frame("x" = rowMeans(counts_filt), # Calculate Mean & Var/Mean
                           "y" = rowVars(counts_filt) / rowMeans(counts_filt)) 
    g_var_all <- rownames(subset(df, x > mn_cut & y > var_cut)) # Genes to keep
    # Plot variance and mean
    par(mfrow = c(1, 1)) # Set the layout to a single plot in a 1x1 grid
    with(df, plot(x, y, main = paste("All Var vs Mean, ", length(g_var_all), 
                                     " genes w/ Mean > ", mn_cut, " & Var > ", var_cut, sep = ""), 
                  pch = 20, cex = 1, xlab = "Mean across all samples", 
                  ylab = "Var/Mean across All samples"))
    with(subset(df, x > mn_cut & y > var_cut), 
         points(x, y, pch = 20, col = "red3", cex = 1)) # Genes above cutoff are red
    abline(v = mn_cut, col = "red3", lty = 2, lwd = 1.5) # Line for mean cutoff
    abline(h = var_cut, col = "red3", lty = 2, lwd = 1.5) # Line for variance cutoff
    removed <- nrow(counts_filt) - length(g_var_all) # Number of filtered genes
    print(paste(removed, "genes removed due to expression below the variance filter"))
    counts_filt <- counts_filt[g_var_all,]  # Update counts_filt 
    #   if (length(unique(groups)) == 2) { # Additional variance by group 
    #     num_groups <- length(unique(groups)) # Determine the number of groups 
    #     # ++++ This needs to be finished still
    #     
    #     removed <- nrow(counts_filt) - length(g_var_group) # Number of filtered genes 
    #     print(paste(removed, " genes removed by the pairwise group variance filter")) 
    #     counts_filt <- counts_filt[g_var_group,]  # Update counts_filt 
    # }
  }
  return(counts_filt) # Final function 1 output
}


