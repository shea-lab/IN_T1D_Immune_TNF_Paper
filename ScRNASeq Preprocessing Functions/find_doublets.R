#' Run DoubletFinder on a Seurat object 
#' 
#' @param seu_object input seurat object 
#' @param sample_name sample name e.g. NP2DPI
#' @param data_dir path to directory of experiment 
#' @param sample_dir path to specific directory where sample is found e.g. SH-1 
#' @param save_opt set this param to TRUE if you want to save the filtered object 
#' 
#' @return DoubletFinder adjusted Seurat object


find_doublets <- function(seu_object, sample_name, data_dir, sample_dir, save_opt = F) {
  # load libraries 
  library(Seurat)
  library(parallel)
  library(DoubletFinder)
  
  # prep for parallelization
  set.seed(8)
  options(mc.cores = detectCores() - 1)
  
  
  ### STEP 1: Data processing ### 
  # Pre-process seurat object with standard seurat workflow
  mouse.sample <- NormalizeData(seu_object)
  mouse.sample <- FindVariableFeatures(mouse.sample)
  mouse.sample <- ScaleData(mouse.sample)
  mouse.sample <- RunPCA(mouse.sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- mouse.sample[["pca"]]@stdev
  sum.stdv <- sum(mouse.sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  mouse.sample <- RunUMAP(mouse.sample, dims = 1:min.pc)
  mouse.sample <- FindNeighbors(object = mouse.sample, dims = 1:min.pc)              
  mouse.sample <- FindClusters(object = mouse.sample, resolution = 0.1)
  
  
  ### STEP 2: identify pK ### 
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(mouse.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  
  ### STEP 3: generate artificial doublets ### 
  ## Homotypic doublet proportion estimate
  annotations <- mouse.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(mouse.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  
  ### STEP 4: run DoubletFinder ### 
  mouse.sample <- doubletFinder(seu = mouse.sample, 
                                PCs = 1:min.pc, 
                                pK = optimal.pk,
                                nExp = nExp.poi.adj)
  metadata <- mouse.sample@meta.data
  colnames(metadata)[grep("DF.classification", colnames(metadata))] <- "doublet_finder"
  mouse.sample@meta.data <- metadata 
  
  
  
  # subset and save
  mouse.singlets <- subset(mouse.sample, doublet_finder == "Singlet")
  
  # compare before to after 
  seu_object 
  mouse.singlets
  
  # clean environment 
  remove(sweep.stats, sweep.list, mouse.sample,  metadata, bcmvn.max, bcmvn)
  #mouse.filtered removed from above
  if (save_opt == TRUE) {
    saveRDS(mouse.singlets, file = paste0("./Data/remove_doublet_filt_", sample, ".rds"))
  }
  
  return(mouse.singlets)
  
}