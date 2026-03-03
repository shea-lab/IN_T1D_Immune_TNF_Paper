#' Run ddqcR on a Seurat object 
#' 
#' @param sample sample name e.g. NP2DPI
#' @param data_dir path to directory of experiment 
#' @param sample_dir path to specific directory where sample is found e.g. SH-1 
#' @param save_opt set this param to TRUE if you want to save the filtered object 
#' 
#' @return ddqcR adjusted Seurat object


filter_ddqcr <- function(data, sample, data_dir, sample_dir, save_opt = F) {
  # load libraries 
  library(Seurat)
  library(ddqcR)
  
  # run initial QC:
  # remove obvious low quality cells- n_genes < 100, percent_mito > 80
  # clustered cells with res 1.3 
  data <- initialQC(data)
  
  df.qc <- ddqc.metrics(data)
  
  df.qc
  
  data <- filterData(data, df.qc)
  data
  
  if (save_opt == TRUE) {
    saveRDS(data, file = paste0("./Data/ddqcR_filt_", sample, ".rds"))
  }
  
  return(data)
  
}