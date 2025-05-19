
#' Control probe PCA
#'  
#' generates PC dataframe from positive control probes to adjust for technical variation
#'  
#' @param ctrl matrix with samples as columns and positive control probe factors as rows 
#' @param samplesheet dataframe with samples for which you want to calculate control probe PCs. Each row is one sample.
#' The sample names must be in a column called "Basename" and they must overlap with the column names in the ctrl matrix
#'
#' @return dataframe with scaled PCs from positive control probes
#' @export
#'
ctrl_probe_PCA <- function(ctrl, samplesheet){
  ctrl <- ctrl[, colnames(ctrl) %in% samplesheet$Basename]
  ctrl_pca <- prcomp(na.omit(t(ctrl))) # run pca
  ctrlprobe_PCAscores <- as.data.frame(ctrl_pca$x) #extract PCA scores
  ctrlprobe_PCAscores <- as.data.frame(scale(ctrlprobe_PCAscores))
  return(ctrlprobe_PCAscores)
}



#' Cell type PCA
#'
#' @param samplesheet samplesheet with samples as rows, has to contain the cell type proportion information
#' @param cell_prop_cols columns of the cell type proportions in samplesheet
#'
#' @return dataframe with scaled PCs from cell type proportions
#' @export
#'
cell_type_PCA <- function(samplesheet, cell_prop_cols){
  
  cell_type_prop_hepi_red <- samplesheet[, cell_prop_cols]
  cellcounts_pca <- prcomp(cell_type_prop_hepi_red) # run pca
  cellcounts_PCAscores <- as.data.frame(cellcounts_pca$x) #extract PCA scores
  cellcounts_PCAscores <- as.data.frame(scale(cellcounts_PCAscores))
  
  return(cellcounts_PCAscores)
}



#' Ancestry PCA
#'
#' @param SNP_info SNP info (can be obtained from DNA methylation data), samples as columns, SNPs/cgs as rows
#' @param samplesheet dataframe with samples for which you want to calculate the PCs. Each row is one sample.
#' The sample names must be in a column called "Basename" and they must overlap with the column names in SNP_info
#'
#' @return dataframe with scaled PCs from ancestry information (SNPs/SNP info from DNA methylation data)
#' @export
#'
ancestry_PCA <- function(SNP_info, samplesheet){
  snps <- SNP_info[, colnames(SNP_info) %in% samplesheet$Basename]
  snps_pca <- prcomp(na.omit(t(snps))) # run pca
  snpsprobe_PCAscores <- as.data.frame(snps_pca$x) #extract PCA scores
  snpsprobe_PCAscores <- as.data.frame(scale(snpsprobe_PCAscores))
  
  return(snpsprobe_PCAscores)
  
}
