

set.seed(123)

# Libraries
library(readxl)
library(writexl)

# Import
our_DMPs <- as.data.frame(read_xlsx("OneDrive - University of Bergen/PhD/24_Case_Control_project/Export_SAFE/output_F/Tables/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_signif_BACON_BH_extended.xlsx"))
load("OneDrive - University of Bergen/PhD/24_Case_Control_project/annotation/UPDATED_ANNOTATION.RData")
Strom <- read.table("OneDrive - University of Bergen/PhD/14_Downloads/SummaryStatistics/OCD_Strom/daner_OCD_full_wo23andMe_190522", header = TRUE)

our_DMPs$GWAS_lowestP <- NA
our_DMPs$GWAS_highest_logOR <- NA

for (i in 1:nrow(our_DMPs)){
  chr <- as.numeric(gsub("chr", "", our_DMPs$CHR[i]))
  pos <- our_DMPs$MAPINFO[i]
  
  Strom_red <- Strom[which(Strom$CHR == chr & (Strom$BP > pos - 1000) & (Strom$BP < pos + 1000) ), ]
  Strom_red$abslogOR <- abs(log(Strom_red$OR))
  
  if (nrow(Strom_red) > 0){
    our_DMPs$GWAS_lowestP[i] <- min(Strom_red$P)
    our_DMPs$GWAS_highest_logOR[i] <- max(Strom_red$abslogOR)
  } else {
    our_DMPs$GWAS_lowestP[i] <- NA
  }
}


anno_upd <- anno_upd[which(anno_upd$CHR != "chrY"), ]

our_DMPs$p_GWAS_permutation <- NA
our_DMPs$z_GWAS_permutation <- NA



our_DMPs$p_GWAS_permutation <- NA
our_DMPs$z_GWAS_permutation <- NA

for(i in 1:nrow(our_DMPs)){
  print(i)
  cgs_minus_ours <- setdiff(anno_upd$Name, our_DMPs$cgID[i])
  
  random_cpgs <- sample(cgs_minus_ours, 1500)
  
  p_random_cpg <- c()
  
  for (cpg in random_cpgs){
    
    chr <- as.numeric(gsub("chr", "", anno_upd$CHR[which(anno_upd$Name == cpg)]))
    pos <- anno_upd$MAPINFO[which(anno_upd$Name == cpg)]
    
    Strom_red <- Strom[which(Strom$CHR == chr & (Strom$BP > pos - 1000) & (Strom$BP < pos + 1000) ), ]
    
    if (nrow(Strom_red) > 0){
      p_random_cpg <- c(p_random_cpg, min(Strom_red$P))
    }
  }
  
  real_p <- as.numeric(our_DMPs$GWAS_lowestP[i])
  random_p <- p_random_cpg[!is.na(p_random_cpg)]
  
  p_val_perm <- (sum(random_p >= real_p) + 1) / (length(random_p) + 1)
  z_score <- (real_p - mean(random_p))/sd(random_p)
  
  our_DMPs$p_GWAS_permutation[i] <- p_val_perm
  our_DMPs$z_GWAS_permutation[i] <- z_score
  
}


write_xlsx(our_DMPs, path = "OneDrive - University of Bergen/PhD/24_Case_Control_project/Comparison_GWAS//Permutation_test/DMPs_females.xlsx")


