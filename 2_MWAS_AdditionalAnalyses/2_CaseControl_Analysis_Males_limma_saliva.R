
################################################################################
# SETUP
################################################################################


# LOAD PATHS
source("S:/Project/WP-epigenetics/04_Pipeline_PhaseI/2_Paths.R")

# SET WD
setwd("S:/Project/WP-epigenetics/10_CaseCtrl_EWAS/")

# LAMBDA FUNCTION
calculate_lambda <- function(pvalues){
  chisq <- qchisq(1-pvalues, 1)
  lambda <- median(chisq)/qchisq(0.5, 1)
  return(lambda)
}



################################################################################
# LOAD LIBRARIES & MAKE FOLDERS TO EXPORT
################################################################################

suppressMessages({
  
  library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
  library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
  library(shiny)
  
  # LOAD LIBRARIES
  library(readxl)
  library(qqman)
  library(ggplot2)
  library(parallel)
  library(doParallel)
  library(writexl)
  library(car)
  library(sva)
  library(limma)
})


# MAKE FOLDERS
if (!dir.exists("output_M")) dir.create("output_M")
if (!dir.exists("output_M/Figures")) dir.create("output_M/Figures")
if (!dir.exists("output_M/Figures/QQplots")) dir.create("output_M/Figures/QQplots")
if (!dir.exists("output_M/Tables")) dir.create("output_M/Tables")
if (!dir.exists("output_M/RData")) dir.create("output_M/RData")


################################################################################
# IMPORT
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

# LOAD SAMPLESHEET & EUROPEAN INFO
samplesheet <- as.data.frame(read_excel(paste0(dir_gen,"/output/Tables/samplesheet_after_QC.xlsx")))

# LOAD M VALUES & BETA VALUES
load(paste0(dir_gen, "/output/RData/Mvalues_final_combined_XY_males.RData"))
load(paste0(dir_gen, "/output/RData/beta_final_combined_XY_males.RData"))

# CTRL PROBES
load(paste0(dir_gen,"/output/RData/Positive_ctrlprobe_intensities.RData")) #ctrl probes

# CELL TYPE PROPORTIONS
load(paste0(dir_gen,"/output/RData/Cell_type_proportions.RData"))

# SNPs FOR ANCESTRY PCs
load(paste0(dir_gen,"/output/RData/comb_SNPs.RData"))

# GENOTYPING DATA
if (available_genotyping_data == "yes"){
  load(genotyping_path)
}

################################################################################
# GENERAL
################################################################################

# SELECT SAMPLES
samples <- samplesheet[which(samplesheet$Sex == "M" & is.na(samplesheet$Transgender)), ]
if (multiple_samples == "yes"){
  samples <- samples[which(samples$Baseline_sample == "Y"),]
}
sample_names <- intersect(samples$Basename, colnames(Mvalues))
samples <- samples[which(samples$Basename %in% sample_names), ]
print(paste0("Nr of Males included: ", length(samples$Basename))) #226
write_xlsx(samples, path = "output_M/Tables/Final_samplesheet_CaseCtrl.xlsx")

Mvalues <- Mvalues[, samples$Basename]
beta <- beta[, samples$Basename]

# VARIABLES
Age <- scale(as.numeric(samples$Age))[, 1]
smoking <- scale(Mvalues["cg05575921", ])[, 1]
CaseCtrl <- factor(samples$Case_Control, levels = c("Control", "Case"))

# CELL TYPE PCs
cellcounts_df <- cell_type_prop_hepi[, c(1:2, 10)]
cellcounts_df = cellcounts_df[rownames(cellcounts_df) %in% sample_names, ]

cell_PCs = prcomp(cellcounts_df)
cell_PCs = as.data.frame(cell_PCs$x)
colnames(cell_PCs) <- paste0("cell_", colnames(cell_PCs))
cell_PCs <- sapply(cell_PCs, scale)



# CTRL PROBE PCs
ctrl <- ctrl[, colnames(ctrl) %in% sample_names]
pca <- prcomp(na.omit(t(ctrl))) # run pca
ctrlprobe_PCAscores <- as.data.frame(pca$x) #extract PCA scores
colnames(ctrlprobe_PCAscores) <- paste0("Ctrl_", colnames(ctrlprobe_PCAscores))
ctrlprobe_PCAscores <- sapply(ctrlprobe_PCAscores, scale)

#### ANCESTRY PCs ###
# PCA

if (available_genotyping_data == "yes"){
  
  if (!sample_names %in% colnames(genotype_matrix)){
    stop("PIPELINE STOPPED. Not all sample names are available in colnames(genotype_matrix). Please either ensure overlapping sample names (Basename in samplesheet) or use available_genotyping_data == no in 2_Paths.R")
  }
  
  genotype_matrix_red <- genotype_matrix[, colnames(genotype_matrix) %in% sample_names]
  pc <- prcomp(t(genotype_matrix_red))
}

if (available_genotyping_data == "no"){
  comb_SNPs_red <- comb_SNPs[, colnames(comb_SNPs) %in% sample_names]
  pc <- prcomp(t(comb_SNPs_red))
}


anc_PCs <- as.data.frame(pc$x)
colnames(anc_PCs) <- paste0("anc_", colnames(anc_PCs))
anc_PCs <- sapply(anc_PCs, scale)

# MAKE DF WITH VARIABLES FOR MODEL
variables_df <- data.frame(CaseCtrl = CaseCtrl, Age = Age, smoking = smoking)
variables_df <- cbind(variables_df, cell_PCs[,1:2], ctrlprobe_PCAscores[,1:15], anc_PCs[,1:10])               


# TEST MULTI-COLINEARITY
temp_variables_df <- variables_df
temp_variables_df$Mvalues <- Mvalues[1, ]

vif_obj <- as.data.frame(vif(lm(Mvalues ~ CaseCtrl + Age + smoking + 
                                  cell_PC1 + cell_PC2 + 
                                  Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 + 
                                  Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                  Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                                  anc_PC1 + anc_PC2 + anc_PC3 + anc_PC4 + anc_PC5 + anc_PC6 + 
                                  anc_PC7 + anc_PC8 + anc_PC9 + anc_PC10, data = temp_variables_df)))
vif_obj$Variable <- rownames(vif_obj)
colnames(vif_obj) <- c("vif", "variable")
write_xlsx(vif_obj, path = "output_M/Tables/VIF_CaseCtrl_Males.xlsx")
rm(temp_variables_df)

high_vif <- vif_obj[6:20, ]
high_vif <- high_vif$variable[which(high_vif$vif > 10)]

################################################################################
# RUN MODEL
################################################################################

# SELECT M VALUES
Mvalues <- Mvalues[row.names(Mvalues) != "cg05575921", ] # remove cpg that is used to adjust for smoking

# TABLE TO SAVE ALL LAMBDA VALUES
lambda_table <- data.frame(Model = c("2Ctrl_2Anc", "5Ctrl_2Anc", "10Ctrl_2Anc", "15Ctrl_2Anc",
                                     "2Ctrl_5Anc", "5Ctrl_5Anc", "10Ctrl_5Anc", "15Ctrl_5Anc",
                                     "2Ctrl_10Anc", "5Ctrl_10Anc", "10Ctrl_10Anc", "15Ctrl_10Anc"), lambda = NA)


### 2 Ctrl PCs, 2 Anc PCs ###
model_variables <- variables_df[, c(1:7, 21:22)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_2CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_2CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_2CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[1,2] <- calculate_lambda(DMPs$P.Value)



### 5 Ctrl PCs, 2 Anc PCs ###
model_variables <- variables_df[, c(1:10, 21:22)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_5CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_5CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_5CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[2,2] <- calculate_lambda(DMPs$P.Value)




### 10 Ctrl PCs, 2 Anc PCs ###
model_variables <- variables_df[, c(1:15, 21:22)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_10CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_10CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_10CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[3,2] <- calculate_lambda(DMPs$P.Value)



### 15 Ctrl PCs, 2 Anc PCs ###
model_variables <- variables_df[, c(1:20, 21:22)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_15CtrlPCs_2AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_15CtrlPCs_2AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_15CtrlPCs_2AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[4,2] <- calculate_lambda(DMPs$P.Value)









### 2 Ctrl PCs, 5 Anc PCs ###
model_variables <- variables_df[, c(1:7, 21:25)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_2CtrlPCs_5AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_2CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_2CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[5,2] <- calculate_lambda(DMPs$P.Value)



### 5 Ctrl PCs, 5 Anc PCs ###
model_variables <- variables_df[, c(1:10, 21:25)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_5CtrlPCs_5AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_5CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_5CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[6,2] <- calculate_lambda(DMPs$P.Value)




### 10 Ctrl PCs, 5 Anc PCs ###
model_variables <- variables_df[, c(1:15, 21:25)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_10CtrlPCs_5AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_10CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_10CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[7,2] <- calculate_lambda(DMPs$P.Value)



### 15 Ctrl PCs, 5 Anc PCs ###
model_variables <- variables_df[, c(1:20, 21:25)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_15CtrlPCs_5AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_15CtrlPCs_5AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_15CtrlPCs_5AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[8,2] <- calculate_lambda(DMPs$P.Value)










### 2 Ctrl PCs, 10 Anc PCs ###
model_variables <- variables_df[, c(1:7, 21:30)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_2CtrlPCs_10AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_2CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_2CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[9,2] <- calculate_lambda(DMPs$P.Value)



### 5 Ctrl PCs, 10 Anc PCs ###
model_variables <- variables_df[, c(1:10, 21:30)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_5CtrlPCs_10AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_5CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_5CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[10,2] <- calculate_lambda(DMPs$P.Value)




### 10 Ctrl PCs, 10 Anc PCs ###
model_variables <- variables_df[, c(1:15, 21:30)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_10CtrlPCs_10AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_10CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_10CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[11,2] <- calculate_lambda(DMPs$P.Value)



### 15 Ctrl PCs, 10 Anc PCs ###
model_variables <- variables_df[, c(1:20, 21:30)]
model_variables <- model_variables[, !names(model_variables) %in% high_vif]

model_formula <- as.formula(paste0(" ~ ", paste0(colnames(model_variables), collapse = " + ")))
print(model_formula)

design <- model.matrix(model_formula, data = variables_df)
fit <- lmFit(Mvalues, design)
fit <- eBayes(fit)

DMPs <- topTable(fit, num=Inf, coef=2, confint = TRUE)
DMPs <- DMPs[order(rownames(DMPs)), ]
DMPs$SE <- (DMPs$CI.R - DMPs$CI.L)/3.92

# EXPORT
#save(DMPs, file = "output_M/RData/DMPs_CaseCtrl_Males_15CtrlPCs_10AncPCs_limma.RData")

pdf(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_15CtrlPCs_10AncPCs_limma.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output_M/Figures/QQplots/QQplot_CaseCtrl_Males_15CtrlPCs_10AncPCs_limma.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

lambda_table[12,2] <- calculate_lambda(DMPs$P.Value)







# EXPORT LAMBDA TABLE
write_xlsx(lambda_table, path = "output_M/Tables/Lambda_CaseCtrl_Males_limma.xlsx")

