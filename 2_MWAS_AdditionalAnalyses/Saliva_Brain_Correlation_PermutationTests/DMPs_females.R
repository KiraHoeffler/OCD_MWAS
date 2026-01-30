

set.seed(123)

# Libraries
library(readxl)
library(writexl)

# Import
our_DMPs <- as.data.frame(read_xlsx("Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/Export_SAFE/output_F/Tables/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_signif_BACON_BH_extended.xlsx"))
corr_table <- read.table("/Users/kirahoeffler/Library/CloudStorage/OneDrive-UniversityofBergen/PhD/14_Downloads/ImageCpG/Illimina_EPIC_covar", header = TRUE, sep = " ", stringsAsFactors = FALSE)
load("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/annotation/UPDATED_ANNOTATION.RData")

extr_our_DMPs <- corr_table[which(corr_table$cgid %in% our_DMPs$cgID), ]

our_DMPs$brain_rho <- NA
our_DMPs$brain_p <- NA
our_DMPs$brain_sign_any <- NA
for (i in 1:nrow(our_DMPs)){
  cpgs <- our_DMPs$cgID[i]
  extr_cgs <- corr_table[which(corr_table$cgid %in% cpgs), ]
  sign_any <- any(extr_cgs$p_rho_brain_saliva <= 0.05)
  
  our_DMPs$brain_rho[i] <- paste0(round(extr_cgs$rho_brain_saliva, 3), collapse = ";")
  our_DMPs$brain_p[i] <- paste0(round(extr_cgs$p_rho_brain_saliva, 3), collapse = ";")
  our_DMPs$brain_sign_any[i] <- sign_any
}

our_DMPs$p_saliva_brain_permutation <- NA
our_DMPs$z_saliva_brain_permutation <- NA

for(i in 1:nrow(our_DMPs)){
  print(i)
  cgs_minus_ours <- setdiff(anno_upd$Name, our_DMPs$cgID[i])
  
  random_cpgs <- sample(cgs_minus_ours, 1500)
  
  real_corr <- abs(as.numeric(our_DMPs$brain_rho[i]))
  random_corr <- abs(corr_table$rho_brain_saliva[which(corr_table$cgid %in% random_cpgs)])
  random_corr <- random_corr[!is.na(random_corr)]
  
  
  p_val_perm <- (sum(random_corr >= real_corr) + 1) / (length(random_corr) + 1)
  z_score <- (real_corr - mean(random_corr))/sd(random_corr)
  
  our_DMPs$p_saliva_brain_permutation[i] <- p_val_perm
  our_DMPs$z_saliva_brain_permutation[i] <- z_score
}

write_xlsx(our_DMPs, path = "Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/Saliva_brain/Permutation_test/DMPs_females.xlsx")
