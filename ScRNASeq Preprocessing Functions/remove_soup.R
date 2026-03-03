#' Run SoupX using CellRanger Outputs
#' 
#' @param sample sample name e.g. NP2DPI
#' @param data_dir path to directory of experiment 
#' @param sample_dir path to specific directory where sample is found e.g. SH-1 
#' 
#' @return SoupX adjusted Seurat object
#' @importFrom Seurat Read10X


remove_soup <- function(sample, data_dir, sample_dir) {
  # load libraries 
  library(Seurat)
  library(SoupX)
  
  # 
  toc = Seurat::Read10X(data.dir = paste0(data_dir, sample_dir, "/sample_filtered_feature_bc_matrix"))
  tod = Seurat::Read10X(data.dir = paste0(data_dir, sample_dir, "/sample_raw_feature_bc_matrix"))
  # tod = Seurat::Read10X("/Users/kategriffin/Documents/Shea Lab/NP Mechanisms in SCI/Flex_Seq/10340-SH/10340-SH-1/count/sample_raw_feature_bc_matrix")
  
  good_cells <-c(toc@Dimnames[[1]])
  tod_good <- tod[good_cells,]
  
  sc = SoupChannel(tod_good, toc)
  
  # Create Seurat object for clustering information
  srat <- CreateSeuratObject(counts = toc)
  srat
  
  srat <- SCTransform(srat, verbose = F)
  srat    <- RunPCA(srat, verbose = F)
  srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat    <- FindClusters(srat, verbose = T)
  
  
  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings
  sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc  <- setDR(sc, umap)
  head(meta)
  
  # auto estimate the soup
  sc  <- autoEstCont(sc)
  
  # output adjusted matrix 
  adj.matrix  <- adjustCounts(sc, roundToInt = T)
  
  # if (file.exists(paste0("./Data/soupX_", sample, "_filt"))) {
  #   user_input <- readline(prompt = "Do you want to rewrite the file? y / n")
  # }
  #   if (user_input == "y" | user_input == "y ") {
  #     file.remove(paste0("./Data/soupX_", sample, "_filt"))
  #     DropletUtils:::write10xCounts(paste0("./Data/soupX_", sample, "_filt"), adj.matrix)
  #   }
  # else {
  #   # save adjusted matrix to use
  #   DropletUtils:::write10xCounts(paste0("./Data/soupX_", sample, "_filt"), adj.matrix)
  # }
  
  
  # output matrix as Seurat object 
  out = CreateSeuratObject(adj.matrix, min.cells = 3)
  
  rm(adj.matrix, meta, sc, srat, toc, tod, tod_good, umap)
  
  return(out)
}
