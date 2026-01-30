
set.seed(123)

# Libraries
library(regioneR)
library(readxl)
library(writexl)

# Import
load("OneDrive - University of Bergen/PhD/24_Case_Control_project/Comparison_GWAS/Non_DMR_regions.RData")
our_DMRs <- as.corr_table.frame(read_xlsx("OneDrive - University of Bergen/PhD/24_Case_Control_project/Export_SAFE/output_meta/Tables/DMRs_CaseCtrl_extended.xlsx"))
corr_table <- read.table("/Users/kirahoeffler/Library/CloudStorage/OneDrive-UniversityofBergen/PhD/14_Downloads/ImageCpG/Illimina_EPIC_covar", header = TRUE, sep = " ", stringsAsFactors = FALSE)
load("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/annotation/UPDATED_ANNOTATION.RData")

genome <- getGenomeAndMask(data.frame(CHR = paste0("chr", c(1:22, "X", "Y")), 
                         chr_length = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                         159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                         115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                         59128983, 63025520, 48129895, 51304566, 155270560, 59373566)))

our_DMRs$brain_rho <- NA
our_DMRs$brain_p <- NA
our_DMRs$brain_sign_any <- NA
our_DMRs$brain_median_corr <- NA
for (i in 1:nrow(our_DMRs)){
  cpgs <- unlist(strsplit(our_DMRs$probe[i], split = ";"))
  extr_cgs <- corr_table[which(corr_table$cgid %in% cpgs), ]
  sign_any <- any(extr_cgs$p_rho_brain_saliva <= 0.05)
  
  our_DMRs$brain_rho[i] <- paste0(round(extr_cgs$rho_brain_saliva, 3), collapse = ";")
  our_DMRs$brain_p[i] <- paste0(round(extr_cgs$p_rho_brain_saliva, 3), collapse = ";")
  our_DMRs$brain_sign_any[i] <- sign_any
  our_DMRs$brain_median_corr[i] <- abs(median(extr_cgs$rho_brain_saliva, na.rm = TRUE))
}
                         
                         
our_DMRs$p_saliva_brain_permutation <- NA
our_DMRs$z_saliva_brain_permutation <- NA

for (j in 1:nrow(our_DMRs)){
  print(j)
  # Convert to GRanges
  DMRs_gr <- toGRanges(our_DMRs[j, c("chr", "start", "end")])
  noDMR_gr <- toGRanges(non_dmrs[, c("CHR", "start", "end")])
  
  random_list <- replicate(
    1500,
    randomizeRegions(
      A = DMRs_gr,
      genome = genome$genome, 
      mask = noDMR_gr
    ),
    simplify = FALSE
  )
  
  rand_median <- c()
  for (i in 1:length(random_list)){
    regions <- as.data.frame(random_list[[i]], row.names = NULL)
    
    cpgs <- anno_upd$Name[which(anno_upd$CHR == as.character(regions$seqnames[1]) & 
                                  anno_upd$MAPINFO >= regions$start[1] & 
                                  anno_upd$MAPINFO <= regions$end[1])]
    extr_cgs <- corr_table[which(corr_table$cgid %in% cpgs), ]
    brain_median_corr <- abs(median(extr_cgs$rho_brain_saliva, na.rm = TRUE))
    rand_median <- c(rand_median, brain_median_corr)
  }
  
  
  real_median <- our_DMRs$brain_median_corr[j]
  rand_median <- rand_median[!is.na(rand_median)]
  
  p_val_perm <- (sum(rand_median >= real_median)) / (length(rand_median)) #0.0579096
  z_score <- (real_median - mean(rand_median))/sd(rand_median)
  
  our_DMRs$p_saliva_brain_permutation[j] <- p_val_perm
  our_DMRs$z_saliva_brain_permutation[j] <- z_score
  
}

write_xlsx(our_DMRs, path = "OneDrive - University of Bergen/PhD/24_Case_Control_project/Saliva_brain/Permutation_test/DMRs_meta.xlsx")
