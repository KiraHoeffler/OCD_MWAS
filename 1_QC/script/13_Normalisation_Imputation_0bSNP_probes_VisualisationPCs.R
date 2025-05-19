

################################################################################
# SETUP
################################################################################

print("10 - NORMALISATION & IMPUTATION - 0bpSNP PROBES")
print("Setting up and loading libraries")

# LOAD PATHS
source("2_Paths.R")

################################################################################
# LOAD LIBERIES
################################################################################

suppressMessages({
  
  #library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
  #library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
  #library(shiny)
  
  # LOAD LIBRARIES
  library(minfi)
  library(limma)
  library(ENmix)
  library(wateRmelon)
  library(ggplot2)
  library(ChAMP)
  library(readxl)
  library(writexl)
  library(factoextra)
  library(plotly)
  library(htmlwidgets)
  
  
  if (array.type == "EPICv2"){
    array <- "IlluminaHumanMethylationEPICv2" #array
    anno_type <- "20a1.hg38" #select type of annotation
    library(IlluminaHumanMethylationEPICv2manifest)
    library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  }
  
  if (array.type == "EPICv1"){
    array <- "IlluminaHumanMethylationEPIC" #array
    anno_type <- "ilm10b4.hg19" #select type of annotation
    library(IlluminaHumanMethylationEPICmanifest)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }
  if (array.type == "450K"){
    array <- "IlluminaHumanMethylation450k" #array
    anno_type <- "ilmn12.hg19" #select type of annotation
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }
})


################################################################################
# IMPORT & ADAPT
################################################################################

# SET PATHS
setwd(dir_gen)

# IMPORT SAMPLESHEET
samplesheet <- as.data.frame(read_excel("output/Tables/Samplesheet.xlsx"))

# LOAD INTENSITIES & DETECTION P VALUES
load('output/RData/intensities_filtered_1e16_SNP0bp.RData')

# LOAD RS SNP PROBES BETA VALUES
load("output/RData/SNP_betas.RData")

# CTRL PROBES
load("output/RData/Positive_ctrlprobe_intensities.RData") #ctrl probes

# CELL TYPE PROPORTIONS
if (tissue_type == "blood"){
  load("output/RData/Cellproportions.RData") 
}
if (tissue_type == "saliva"){
  load("output/RData/Cell_type_proportions.RData")
}


# THEMES
load("output/RData/theme.RData")
load("output/RData/theme_transparent.RData")

# FUNCTION: calculate beta from genotypes (from meffil package)
calculate.beta.genotypes <- function(snp.betas, centers=c(0.2,0.5,0.8)) {
  x <- t(apply(snp.betas,1,function(x) {
    tryCatch(kmeans(x, centers=centers)$cluster - 1,
             error=function(e) {
               cluster <- rep(1,ncol(snp.betas))
               cluster[which(x < min(centers))] <- 0
               cluster[which(x > max(centers))] <- 2
               cluster
             })
  }))
  dimnames(x) <- dimnames(snp.betas)
  x
}

################################################################################
# NORMALISE IN INTENSITIES & CALCULATE BETA VALUE
################################################################################

print("Step 1/4: Quantile normalisation of intensities")

# NORMALISE
TypeII.Green = normalizeQuantiles(TypeII.Green)
TypeII.Red = normalizeQuantiles(TypeII.Red)
TypeI.Green.M = normalizeQuantiles(TypeI.Green.M)
TypeI.Green.U = normalizeQuantiles(TypeI.Green.U)
TypeI.Red.M = normalizeQuantiles(TypeI.Red.M)
TypeI.Red.U = normalizeQuantiles(TypeI.Red.U)
save(TypeII.Green, TypeII.Red, TypeI.Green.M, TypeI.Green.U, TypeI.Red.M, TypeI.Red.U, file = "output/RData/intensities_normalised_SNP0bp.RData")

# CALCULATE BETA VALUES & EXPORT
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta_SNP0bp = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))

rm(TypeII.Green,TypeII.Red,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U,TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)


################################################################################
# IMPUTE MISSING VALUES
################################################################################

print("Step 2/4: impute missing values.")

samplesheet <- samplesheet[order(samplesheet$Basename), ]
samplesheet <- samplesheet[which(samplesheet$Basename %in% colnames(beta_SNP0bp)), ]

suppressMessages({
  invisible(capture.output({
    imputed_beta <- champ.impute(beta = beta_SNP0bp, pd = samplesheet, method = "KNN")
  }))
})

beta_SNP0bp_imp <- imputed_beta$beta

save(beta_SNP0bp_imp, file='output/RData/beta_QN_1e16_SNP0bp_KNN.RData',compress=F)


# VISUALISE BETA DISTRIBUTION
pdf("output/Figures/Distributions/Beta_distribution_filtered_QN_KNN_SNP0bp.pdf", width = 7, height = 5)
densityPlot(beta_SNP0bp_imp, main = "Normalised & imputed SNP0bp beta values")
dev.off()



################################################################################
# CALL GENOTYPES OF SNP PROBES
################################################################################

print("Step 3/4: Call genotypes of rs SNP probes.")

# exclude samples from rs SNPs that were excluded in filtering step
snp_betas <- snp_betas[, colnames(snp_betas) %in% colnames(beta_SNP0bp_imp)]

rs_genotypes <- calculate.beta.genotypes(snp_betas)
rs_genotypes[rs_genotypes == 1] <- 0.5
rs_genotypes[rs_genotypes == 2] <- 1

save(rs_genotypes, file = "output/RData/rs_genotypes.RData")



################################################################################
# VISUALISE CELL TYPE PCs, CONTROL PROBE PCs, AND ANCESTRY PCs
################################################################################

print("Step 4/4: Visualise cell type PCs, control probe PCs, and ancestry PCs.")

### EXTRACT FROM SAMPLESHEET ###
Age <- samplesheet$Age
Sex <- samplesheet$Sex


### CONTROL PROBE PCA ###
ctrl <- ctrl[, colnames(ctrl) %in% samplesheet$Basename]
ctrl_pca <- prcomp(na.omit(t(ctrl))) # run pca

ctrlprobe_PCAscores <- as.data.frame(ctrl_pca$x) #extract PCA scores

Ctrl_PC1 <- scale(ctrlprobe_PCAscores$PC1)[, 1]
Ctrl_PC2 <- scale(ctrlprobe_PCAscores$PC2)[, 1]
Ctrl_PC3 <- scale(ctrlprobe_PCAscores$PC3)[, 1]
Ctrl_PC4 <- scale(ctrlprobe_PCAscores$PC4)[, 1]
Ctrl_PC5 <- scale(ctrlprobe_PCAscores$PC5)[, 1]
Ctrl_PC6 <- scale(ctrlprobe_PCAscores$PC6)[, 1]
Ctrl_PC7 <- scale(ctrlprobe_PCAscores$PC7)[, 1]
Ctrl_PC8 <- scale(ctrlprobe_PCAscores$PC8)[, 1]
Ctrl_PC9 <- scale(ctrlprobe_PCAscores$PC9)[, 1]
Ctrl_PC10 <- scale(ctrlprobe_PCAscores$PC10)[, 1]

# Scree plot to see how much variance the PCs explain
scree <- fviz_eig(ctrl_pca, main = "Control probe PCs")
ggsave("output/Figures/Scree_CtrlProbes.svg", scree, device = "svg", units = "cm", width = 10, height = 8, dpi = 400)

ctrl_df <- scree$data
write_xlsx(ctrl_df, path = "output/Tables/Scree_CtrlProbes.xlsx")


### CELL TYPE PCA ###
if (tissue_type == "blood"){
  if (array.type == "450K" | array.type == "EPICv1"){
    cellcounts <- cellcounts$prop
  }
  
  cellcounts <- cellcounts[rownames(cellcounts) %in% samplesheet$Basename, ]
  cellcounts_pca <- prcomp(na.omit(cellcounts)) # run pca
  
  cellcounts_PCAscores <- as.data.frame(cellcounts_pca$x) #extract PCA scores
  
  cellcounts_PC1 <- scale(cellcounts_PCAscores$PC1)[, 1]
  cellcounts_PC2 <- scale(cellcounts_PCAscores$PC2)[, 1]
  cellcounts_PC3 <- scale(cellcounts_PCAscores$PC3)[, 1]
  cellcounts_PC4 <- scale(cellcounts_PCAscores$PC4)[, 1]
  cellcounts_PC5 <- scale(cellcounts_PCAscores$PC5)[, 1]
  
  # Scree plot to see how much variance the PCs explain
  scree <- fviz_eig(cellcounts_pca, main = "Cell count PCs")
  ggsave("output/Figures/Scree_cellcounts.svg", scree, device = "svg", units = "cm", width = 10, height = 8, dpi = 400)
  
  cellcounts_df <- scree$data
  write_xlsx(cellcounts_df, path = "output/Tables/Scree_cellcounts.xlsx")
  
  
  ### RESIDUALISED BETA SNP0bp ###
  beta_SNP0bp_imp <- beta_SNP0bp_imp[, colnames(beta_SNP0bp_imp) %in% samplesheet$Basename]
  resid_beta <- as.data.frame(matrix(data = NA, nrow = nrow(beta_SNP0bp_imp), ncol = ncol(beta_SNP0bp_imp)))
  
  for (i in c(1:nrow(beta_SNP0bp_imp))){
    model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + cellcounts_PC3 + cellcounts_PC4 + cellcounts_PC5 +
                  Sex + Age + 
                  Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
    residuals <- resid(model)
    resid_beta[i, ] <- residuals
  }
  rownames(resid_beta) <- rownames(beta_SNP0bp_imp)
  colnames(resid_beta) <- colnames(beta_SNP0bp_imp)
}

if (tissue_type == "saliva"){
  cellcounts <- data.frame(Epi = samplesheet$Epi, Fib = samplesheet$Fib, comb_ICs = samplesheet$comb_ICs)
  cellcounts_pca <- prcomp(na.omit(cellcounts)) # run pca
  
  cellcounts_PCAscores <- as.data.frame(cellcounts_pca$x) #extract PCA scores
  
  cellcounts_PC1 <- scale(cellcounts_PCAscores$PC1)[, 1]
  cellcounts_PC2 <- scale(cellcounts_PCAscores$PC2)[, 1]
  
  # Scree plot to see how much variance the PCs explain
  scree <- fviz_eig(cellcounts_pca, main = "Cell count PCs")
  ggsave("output/Figures/Scree_cellcounts.svg", scree, device = "svg", units = "cm", width = 10, height = 8, dpi = 400)
  
  cellcounts_df <- scree$data
  write_xlsx(cellcounts_df, path = "output/Tables/Scree_cellcounts.xlsx")
  
  
  ### RESIDUALISED BETA SNP0bp ###
  beta_SNP0bp_imp <- beta_SNP0bp_imp[, colnames(beta_SNP0bp_imp) %in% samplesheet$Basename]
  resid_beta <- as.data.frame(matrix(data = NA, nrow = nrow(beta_SNP0bp_imp), ncol = ncol(beta_SNP0bp_imp)))
  
  for (i in c(1:nrow(beta_SNP0bp_imp))){
    model <- lm(beta_SNP0bp_imp[i, ] ~ cellcounts_PC1 + cellcounts_PC2 + 
                  Sex + Age + 
                  Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10)
    residuals <- resid(model)
    resid_beta[i, ] <- residuals
  }
  rownames(resid_beta) <- rownames(beta_SNP0bp_imp)
  colnames(resid_beta) <- colnames(beta_SNP0bp_imp)
}

### COMBINE RESIDUALISED SNP0bp with rs beta ###

rs_genotypes <- rs_genotypes[, colnames(rs_genotypes) %in% colnames(resid_beta)]

comb_SNPs <- rbind(resid_beta, rs_genotypes)
save(comb_SNPs, file = "output/RData/comb_SNPs.RData")

### PCA ###

pc <- prcomp(t(comb_SNPs))
top10pc <- as.data.frame(pc$x[,1:10])
top10pc$Basename <- rownames(top10pc)

pc_plot <- qplot(x=PC1, y=PC2, data = top10pc, color = samplesheet$Ancestry, ) + th + th_transparent + 
  theme(legend.background = element_rect(color = NA, fill = NA),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank())
ggsave("output/Figures/PCA_Ancestry.svg", pc_plot, units = "cm", device = "svg", width = 13, height = 8, dpi = 400)

# 3D-plot
plot <- plot_ly(top10pc, x = ~PC1, y = ~PC2, z = ~PC3, color = samplesheet$Ancestry,
                type = "scatter3d", mode = "markers", marker = list(size = 4, opacity = 1))
plot <- plot %>% layout(scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3")))
plot

saveWidget(plot, file = "output/Figures/Ancestry_PCs.html")
saveWidget(plot, file = paste0(dir_gen,"Ancestry_PCs.html"))


rm(list = ls());gc()


