
################################################################################
# SETUP
################################################################################

print("12 - COMBINE AUTOSOMES AND X/Y CHR - MALES")
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
})

################################################################################
# IMPORT
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

load("output/RData/beta_final_XY_males.RData")
beta_m <- beta
load("output/RData/beta_final_autosomes.RData")


################################################################################
# COMBINE DF
################################################################################

print("Combine dfs and visualise the distributions.")

beta_m_auto <- beta[, colnames(beta) %in% colnames(beta_m)]
rm(beta)

beta <- rbind(beta_m_auto, beta_m)
print(paste0("Nr of males in finale combined male dataset: ", ncol(beta)))


# SAVE FINAL BETA VALUE DF
save(beta, file = "output/RData/beta_final_combined_XY_males.RData")


# REPLACE ALL ZEROS IN BETA VALUE WITH LOWEST VALUE
replacement_zero <- min(beta[which(beta > 0 & !is.na(beta))]) 
beta_temp <- beta
beta_temp[beta_temp == 0] <- replacement_zero

# CALCULATE M VALUES 
Mvalues <- B2M(beta_temp)
rm(beta_temp)

#EXPORT M VALUES
save(Mvalues, file='output/RData/Mvalues_final_combined_XY_males.RData',compress=F)


# VISUALISE BETA DISTRIBUTION
pdf("output/Figures/Distributions/beta_distribution_final_comb_XY_males.pdf", width = 7, height = 5)
densityPlot(beta, main = "Filtered and normalised beta values")
dev.off()

# VISUALISE M VALUE DISTRIBUTION
pdf("output/Figures/Distributions/Mvalue_distribution_final_comb_XY_males.pdf", width = 7, height = 5)
densityPlot(Mvalues, main = "Filtered and normalised M values")
dev.off()



rm(list = ls());gc()

