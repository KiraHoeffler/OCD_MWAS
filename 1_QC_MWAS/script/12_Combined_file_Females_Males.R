
################################################################################
# SETUP
################################################################################

print("13 - COMBINE AUTOSOMES + X CHR - MALES AND FEMALES")
print("Setting up and loading libraries")

source("2_Paths.R")


################################################################################
# LIBRARIES
################################################################################

suppressMessages({
  #library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
  #library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
  #library(shiny)
  
  library(ENmix)
  library(minfi)
  library(readxl)
  library(writexl)
})


################################################################################
# IMPORT
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

load("output/RData/beta_final_combined_XY_males.RData")
beta_male <- beta
load("output/RData/beta_final_combined_X_females.RData")
beta_female <- beta
rm(beta)

################################################################################
# COMBINE DF
################################################################################

print("Step 1/2: Combine autosomes + X chr for females and males.")

common_cpgs <- intersect(rownames(beta_female), rownames(beta_male))

beta_male <- beta_male[rownames(beta_male) %in% common_cpgs, ]
beta_female <-beta_female[rownames(beta_female) %in% common_cpgs, ]

beta <- cbind(beta_male, beta_female)

# SAVE FINAL BETA VALUE DF
save(beta, file = "output/RData/beta_final_combined_males_females.RData")


# REPLACE ALL ZEROS IN BETA VALUE WITH LOWEST VALUE
replacement_zero <- min(beta[which(beta > 0 & !is.na(beta))]) 
beta_temp <- beta
beta_temp[beta_temp == 0] <- replacement_zero

# CALCULATE M VALUES 
Mvalues <- B2M(beta_temp)
rm(beta_temp)

#EXPORT M VALUES
save(Mvalues, file='output/RData/Mvalues_final_combined_males_females.RData',compress=F)


# VISUALISE BETA DISTRIBUTION
pdf("output/Figures/Distributions/beta_distribution_final_comb_males_females.pdf", width = 7, height = 5)
densityPlot(beta, main = "Filtered and normalised beta values")
dev.off()

# VISUALISE M VALUE DISTRIBUTION
pdf("output/Figures/Distributions/Mvalue_distribution_final_comb_males_females.pdf", width = 7, height = 5)
densityPlot(Mvalues, main = "Filtered and normalised M values")
dev.off()



################################################################################
# REMOVE INDIVIDUALS THAT DID NOT PASS SEX CHROMOSOME QC FROM AUTOSOME DFs
################################################################################

print("Step 2/2: Combine individuals that did not pass sex chr QC from autosome dfs.")

indiv_passing_QC <- colnames(beta)

rm(beta)
rm(Mvalues)
load("output/RData/beta_final_autosomes.RData")
load("output/RData/Mvalues_final_autosomes.RData")

# UPDATE BETA VALUES
beta <- beta[, colnames(beta) %in% indiv_passing_QC]
save(beta, file = "output/RData/beta_final_autosomes.RData")

# UPDATE M VALUES
Mvalues <- Mvalues[, colnames(Mvalues) %in% indiv_passing_QC]
save(Mvalues, file = "output/RData/Mvalues_final_autosomes.RData")

# IMPORT SAMPLESHEET
samplesheet <- as.data.frame(read_excel("output/Tables/Samplesheet.xlsx"))

# SAVE SAMPLESHEET AFTER QC
samples_red <- samplesheet[which(samplesheet$Basename %in% indiv_passing_QC), ]
samples_red <- samples_red[order(samples_red$Basename), ]
write_xlsx(samples_red, path = "output/Tables/samplesheet_after_QC.xlsx")

######
rm(list = ls());gc()
