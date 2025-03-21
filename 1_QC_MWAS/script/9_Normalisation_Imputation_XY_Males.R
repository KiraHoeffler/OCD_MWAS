

################################################################################
# SETUP
################################################################################

print("9 - NORMALISATION & IMPUTATION - X AND Y CHR IN MALES")
print("Setting up and loading libraries")

# LOAD PATHS
source("2_Paths.R")

if (array.type == "EPICv1"){
  annofile_path <- file.path(resources,"EPICv1_manifest.csv") #downloaded Illumina annotation file
}
if (array.type == "EPICv2"){
  annofile_path <- file.path(resources,"EPICv2_manifest.csv")
  replic_path <- file.path(resources,"HandleReplicProbes.RData")
}
if(array.type == "450K"){
  annofile_path <- file.path(resources,"450k_manifest.csv")
}

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

# LOAD INTENSITIES & DETECTION P VALUES
load('output/RData/intensities_filtered_1e16_XY_males.RData')

# ANNOTATION
annotation <- read.csv(annofile_path, skip = 7)

# IMPORT SAMPLESHEET
samplesheet <- as.data.frame(read_excel("output/Tables/Samplesheet.xlsx"))


################################################################################
# QUANTILE-NORMALISATION
################################################################################

print("Step 1/2: Quantile normalisation of intensities")

# NORMALISE
TypeII.Green = normalizeQuantiles(TypeII.Green)
TypeII.Red = normalizeQuantiles(TypeII.Red)
TypeI.Green.M = normalizeQuantiles(TypeI.Green.M)
TypeI.Green.U = normalizeQuantiles(TypeI.Green.U)
TypeI.Red.M = normalizeQuantiles(TypeI.Red.M)
TypeI.Red.U = normalizeQuantiles(TypeI.Red.U)
save(TypeII.Green, TypeII.Red, TypeI.Green.M, TypeI.Green.U, TypeI.Red.M, TypeI.Red.U, file = "output/RData/intensities_normalised_XY_males.RData")

# CALCULATE BETA VALUES & EXPORT
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))

rm(TypeII.Green,TypeII.Red,TypeI.Green.M,TypeI.Green.U,TypeI.Red.M,TypeI.Red.U,TypeII.betas,TypeI.Green.betas,TypeI.Red.betas)
#save(beta, file='output/RData/beta_QN_1e16.RData',compress=F)

# VISUALISE BETA DISTRIBUTION
pdf("output/Figures/Distributions/beta_distribution_filtered_QN_1e16_XY_males.pdf", width = 7, height = 5)
densityPlot(beta, main = "Filtered and QN beta values")
dev.off()




################################################################################
# DISTRIBUTION NA VALUES PER CPG FOR BETA VALUES
################################################################################

print("Step 2/2: impute missing values.")

# COUNT NAs PER CPG
num_na <- as.data.frame(rowSums(is.na(beta)))
colnames(num_na) <- "nr_NA"
write_xlsx(num_na, path = "output/Tables/Nr_NAs_per_CpGs_XY_males.xlsx")

################################################################################
# IMPUTATION OF BETA VALUES
################################################################################

suppressMessages({
  invisible(capture.output({
    imputed_beta <- champ.impute(beta = beta, pd = samplesheet, method = "KNN")
  }))
})

beta_imp <- imputed_beta$beta
sum(is.na(beta_imp)) # should be zero



# EXPORT
if (array.type == "EPICv2"){
  save(beta_imp, file = "output/RData/beta_filtered_QN_BMIQ_imputed_1e16_XY_males.RData", compress=F)
}
if (array.type == "EPICv1" | array.type == "450K"){
  beta <- beta_imp
  save(beta, file = "output/RData/beta_final_XY_males.RData", compress=F)
}


################################################################################
# REPLICATED PROBES (ONLY EPICv2)
################################################################################

print("Extra step EPICv2: Handle replicated probes.")

if (array.type == "EPICv2"){
  load(replic_path)
  
  beta <- as.matrix(HandleReplicProbes(beta_imp))
  
  # SAVE FINAL BETA VALUE DF
  save(beta, file = "output/RData/beta_final_XY_males.RData")
}

################################################################################
# VISUALISE BETA DISTRIBUTION
################################################################################

pdf("output/Figures/Distributions/Beta_distribution_final_XY_males.pdf", width = 7, height = 5)
densityPlot(beta, main = "final beta - XY males")
dev.off()


rm(list = ls());gc()

