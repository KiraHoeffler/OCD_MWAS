
set.seed(123)

# Libraries
library(regioneR)
library(readxl)
library(writexl)

# Import
our_DMRs <- as.data.frame(read_xlsx("OneDrive - University of Bergen/PhD/24_Case_Control_project/Export_SAFE/output_meta/Tables/DMRs_CaseCtrl_extended.xlsx"))
load("OneDrive - University of Bergen/PhD/24_Case_Control_project/annotation/UPDATED_ANNOTATION.RData")
Strom <- read.table("OneDrive - University of Bergen/PhD/14_Downloads/SummaryStatistics/OCD_Strom/daner_OCD_full_wo23andMe_190522", header = TRUE)

our_DMRs$GWAS_lowestP <- NA
our_DMRs$GWAS_highest_logOR <- NA

for (i in 1:nrow(our_DMRs)){
  chr <- as.numeric(gsub("chr", "", our_DMRs$chr[i]))
  start_pos <- our_DMRs$start[i]
  end_pos <- our_DMRs$end[i]
  
  Strom_red <- Strom[which(Strom$CHR == chr & (Strom$BP > start_pos - 1000) & (Strom$BP < end_pos + 1000) ), ]
  Strom_red$abslogOR <- abs(log(Strom_red$OR))
  
  if (nrow(Strom_red) > 0){
    our_DMRs$GWAS_lowestP[i] <- min(Strom_red$P)
    our_DMRs$GWAS_highest_logOR[i] <- max(Strom_red$abslogOR)
  } else {
    our_DMRs$GWAS_lowestP[i] <- NA
  }
}



genome <- getGenomeAndMask(data.frame(CHR = paste0("chr", c(1:22, "X", "Y")), 
                                      chr_length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                                                     159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                                                     115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                                                     59128983, 63025520, 48129895, 51304566, 0, 0)))


our_DMRs$p_GWAS_permutation <- NA
our_DMRs$z_GWAS_permutation <- NA

for (j in 1:nrow(our_DMRs)){
  
  if(!is.na(our_DMRs$GWAS_lowestP[j])){
      print(j)
      # Convert to GRanges
      DMRs_gr <- toGRanges(our_DMRs[j, c("chr", "start", "end")])
      
      random_list <- replicate(
        1500,
        randomizeRegions(
          A = DMRs_gr,
          genome = genome$genome, 
          mask = DMRs_gr
        ),
        simplify = FALSE
      )
      
      rand_pval <- c()
      for (i in 1:length(random_list)){
        regions <- as.data.frame(random_list[[i]], row.names = NULL)
        
        chr <- as.numeric(gsub("chr", "", regions$seqnames[1]))
        start_pos <- regions$start[1]
        end_pos <- regions$end[1]
        
        Strom_red <- Strom[which(Strom$CHR == chr & (Strom$BP > start_pos - 1000) & (Strom$BP < end_pos + 1000) ), ]
        
        if (nrow(Strom_red) > 0){
          rand_pval <- c(rand_pval, min(Strom_red$P))
        } 
      }
      
      
      real_pval <- -log10(our_DMRs$GWAS_lowestP[j])
      rand_pval <- -log10(rand_pval[!is.na(rand_pval)])
      
      p_val_perm <- (sum(rand_pval >= real_pval)) / (length(rand_pval))
      z_score <- (real_pval - mean(rand_pval))/sd(rand_pval)
      
      our_DMRs$p_GWAS_permutation[j] <- p_val_perm
      our_DMRs$z_GWAS_permutation[j] <- z_score
  }
}

write_xlsx(our_DMRs, path = "OneDrive - University of Bergen/PhD/24_Case_Control_project/Comparison_GWAS/Permutation_test/DMR_meta.xlsx")


