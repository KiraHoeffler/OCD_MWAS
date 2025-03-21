
################################################################################
# SETUP
################################################################################

setwd("S:/Project/WP-epigenetics/10_CaseCtrl_EWAS/")

# PATHS
DMPs_path <- "S:/Project/WP-epigenetics/04_Pipeline_PhaseI/output/RData/Sex_adjusted/DMPs_CaseCtrl_Sex_adjusted_15CtrlPCs_2AncPCs_limma.RData"
anno_path <- "S://Project/WP-epigenetics/02_Import/Illumina/EPIC-8v2-0_A1.csv"
beta_path <- "S:/Project/WP-epigenetics/04_Pipeline_PhaseI/output/RData/beta_final_autosomes.RData"
Mvalues_path <- "S:/Project/WP-epigenetics/04_Pipeline_PhaseI/output/RData/Mvalues_final_autosomes.RData"
ctrl_path <- "S://Project/WP-epigenetics/04_Pipeline_PhaseI/output/RData/Positive_ctrlprobe_intensities.RData"
SNP_path <- "S://Project/WP-epigenetics/04_Pipeline_PhaseI/output/RData/comb_SNPs.RData"


# MAKE FOLDERS
if (!dir.exists("output/Figures/DMRs")) dir.create("output/Figures/DMRs")
if (!dir.exists("output/Figures/DMRs/raw")) dir.create("output/Figures/DMRs/raw")
if (!dir.exists("output/Figures/DMRs/combined")) dir.create("output/Figures/DMRs/combined")
if (!dir.exists("output/Figures/DMRs/unmedicated")) dir.create("output/Figures/DMRs/unmedicated")
if (!dir.exists("output/Figures/DMRs/noComorbidity/combined")) dir.create("output/Figures/DMRs/noComorbidity/combined")
if (!dir.exists("output/Figures/DMRs/noComorbidity")) dir.create("output/Figures/DMRs/noComorbidity")
if (!dir.exists("output/Figures/DMRs/unmedicated/combined")) dir.create("output/Figures/DMRs/unmedicated/combined")
if (!dir.exists("output/Figures/DMRs/european")) dir.create("output/Figures/DMRs/european")
if (!dir.exists("output/Figures/DMRs/european/combined")) dir.create("output/Figures/DMRs/european/combined")
if (!dir.exists("output/Figures/DMRs/mQTL")) dir.create("output/Figures/DMRs/mQTL")
if (!dir.exists("output/Figures/DMRs/mQTL/combined")) dir.create("output/Figures/DMRs/mQTL/combined")


################################################################################
# LIBRARIES
################################################################################

library(pkgload, lib.loc = "C:/Users/kihof7027/Downloads/")
library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
library(shiny, lib.loc = "Z:/Bioconductorpackages/319")
library(readxl)
library(bacon)
library(qqman)
library(writexl)
library(dmrff)
library(dplyr)
library(missMethyl)
library(limma)
library(ENmix)
library(patchwork)
library(irlba)
library(bigsnpr)
library(COMBAT)

################################################################################
# LOAD & ADAPT
################################################################################

#samplesheet
samples <- as.data.frame(read_xlsx("output/Tables/Final_samplesheet_CaseCtrl.xlsx"))
mean(samples$Pre_YB_Sum, na.rm = TRUE)
sd(samples$Pre_YB_Sum, na.rm = TRUE)

#DMPS
load(DMPs_path)
DMPs$cgID <- rownames(DMPs)

#functions
source(file = "../05_Analysis_old/0_General/Functions/lambda.R")
source(file = "../05_Analysis_old/0_General/Functions/DMR_functions.R")
source(file = "../05_Analysis_old/0_General/Functions/PCAs_for_model.R")

#annotation
annotation <- read.csv(anno_path, skip = 7)

#themes
load("output/RData/theme.RData")
load("output/RData/theme_transparent.Rdata")

# CTRL PROBES
load(ctrl_path) #ctrl probes

# SNPs FOR ANCESTRY PCs
load(SNP_path)


################################################################################
# ADAPT FINAL MODEL
################################################################################


# APPLY BACON CORRECTION
DMPs_bacon <- bacon::bacon(effectsizes = DMPs$logFC, standarderrors = DMPs$SE)

DMPs <- data.frame(cgID = DMPs$cgID, 
                   effectsize = bacon::es(DMPs_bacon),
                   t = bacon::tstat(DMPs_bacon),
                   P.Value = bacon::pval(DMPs_bacon), 
                   SE = bacon::se(DMPs_bacon))

lambda(DMPs$P.Value)#1.03966


# CORRECT FOR MULTIPLE TESTING
DMPs$adj_p_BH <- p.adjust(DMPs$P.Value, method = "BH")
DMPs$adj_p_Bonferroni <- p.adjust(DMPs$P.Value, method = "bonferroni")

# EXPORT
save(DMPs, file = "output/RData/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_BACON.RData")

pdf(file = "output/Figures/QQplot_CaseCtrl_15CTRLPCs_2AncPCs_limma_BACON.pdf", width = 5, height = 5)
qq(DMPs$P.Value)
dev.off()

png(file = "output/Figures/QQplot_CaseCtrl_15CTRLPCs_2AncPCs_limma_BACON.png", width = 5, height = 5, units = "in", res = 300)
qq(DMPs$P.Value)
dev.off()

# ANNOTATION
anno <- unique(annotation[which(annotation$Name %in% DMPs$cgID), c("Name", "CHR", "MAPINFO")])
anno <- anno[order(anno$Name), ] # order alphabetically

# ANNOTATE RESULTS AND SAVE
DMPs_annotation <- merge(DMPs, annotation, by.x = "cgID", by.y = "Name", all.x = TRUE)
save(DMPs_annotation, file = "output/RData/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_annotated_BACON.RData")
DMPs_anno <- merge(DMPs, anno, by.x = "cgID", by.y = "Name", all.x = TRUE)
save(DMPs_anno, file = "output/RData/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_anno_BACON.RData")
rm(DMPs)

# EXPORT SIGNIFICANT RESULTS
signif_DMPs <- DMPs_annotation[which(DMPs_annotation$adj_p_Bonferroni <= 0.05),]
nrow(signif_DMPs) #0
write_xlsx(signif_DMPs, path = "output/Tables/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_signif_BACON.xlsx")
save(signif_DMPs, file = "output/RData/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_signif_BACON.RData")

signif_DMPs_BH <- DMPs_annotation[which(DMPs_annotation$adj_p_BH <= 0.05),]
nrow(signif_DMPs_BH) #14
write_xlsx(signif_DMPs_BH, path = "output/Tables/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_signif_BACON_BH.xlsx")
save(signif_DMPs_BH, file = "output/RData/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_signif_BACON_BH.RData")


rm(DMPs_annotation)



################################################################################
# DMRs
################################################################################

####################### PREPARE ################################################

# LOAD AND FILTER BETA VALUES
load(beta_path)
beta_filt <- beta[row.names(beta) %in% DMPs_anno$cgID, colnames(beta) %in% samples$Basename]

########################## DMRFF ###############################################

dmrff_result <- dmrff(estimate = DMPs_anno$effectsize,
                      se = DMPs_anno$SE,
                      p.value = DMPs_anno$P.Value,
                      methylation = beta_filt,
                      chr = DMPs_anno$CHR,
                      pos = DMPs_anno$MAPINFO,
                      maxgap = 500,
                      p.cutoff = 0.05)
#41501     candidate regions

save(dmrff_result, file = "output/RData/DMRs_dmrff_CaseCtrl_Females_15CtrlPCs_2AncPCs_BACON.RData")


############################ FILTER AND ANNOTATE DMRFF##########################

nr_signif_CpGs <- 3

# only select significant DMRs and those with at least 3 CpGs
dmrff_sign <- dmrff_result[which(dmrff_result$p.adjust <= 0.05 & dmrff_result$n >= 3), ] 

# nominal significant DMPs
DMPs_anno_sign05 <- DMPs_anno[which(DMPs_anno$P.Value < 0.05), ]

#annotate DMP_overlap
dmrff_sign$GeneString <- NA
dmrff_sign$Nr_CpGs <- NA
dmrff_sign$Nr_signif_CpGs_05 <- NA


for (i in c(1:nrow(dmrff_sign))){
  chr <- dmrff_sign$chr[i]
  start <- dmrff_sign$start[i]
  end <- dmrff_sign$end[i]
  
  GeneStart <- unique(annotation$UCSC_RefGene_Name[which(annotation$CHR == chr & annotation$MAPINFO == start)])
  GeneEnd <- unique(annotation$UCSC_RefGene_Name[which(annotation$CHR == chr & annotation$MAPINFO == end)])
  Genes <- unique(c(unlist(strsplit(GeneStart, ";")), unlist(strsplit(GeneEnd, ";"))))
  GeneString <- paste(Genes, collapse = ";")
  dmrff_sign$GeneString[i] <- GeneString
  
  DMR_cpgs <- unique(annotation$Name[which(annotation$CHR == chr & annotation$MAPINFO >= start & annotation$MAPINFO <= end)])
  dmrff_sign$Nr_CpGs[i] <- length(DMR_cpgs)
  nr_sign05 <- length(intersect(DMR_cpgs, DMPs_anno_sign05$cgID))
  dmrff_sign$Nr_signif_CpGs_05[i] <- nr_sign05
}

dmrff_sign <- dmrff_sign[which(dmrff_sign$Nr_signif_CpGs_05 >= nr_signif_CpGs), ]
dmrff_sign <- dmrff_sign[which(dmrff_sign$p.adjust <= 0.001), ]

dmrff_sign$Perc_sign_CpGs <- dmrff_sign$Nr_signif_CpGs_05 / (dmrff_sign$Nr_CpGs)
dmrff_sign <- dmrff_sign[which(dmrff_sign$Perc_sign_CpGs > 0.5), ]

print(paste0("significant DMRs after filtering: ", nrow(dmrff_sign))) #33

# order: adj. p value
dmrff_sign <- dmrff_sign[order(dmrff_sign$p.adjust), ]

# add UNANNOTATED
length_unannotated <- length(dmrff_sign$GeneString[which(dmrff_sign$GeneString == "")])

count <- 1
for (i in c(1:nrow(dmrff_sign))){
  if (dmrff_sign$GeneString[i] == ""){
    New_Gene <- paste0("UNANNOTATED", count)
    dmrff_sign$GeneString[i] <- New_Gene
    count <- count + 1
  }
}


write_xlsx(dmrff_sign, "output/Tables/DMRs_dmrff_CaseCtrl_15CtrlPCs_2AncPCs_filtered_BACON.xlsx")

# EXTRACT CPGS IN DMRs
all_DMR_cpgs <- extract_CpGs_from_DMRs(dmrff_sign, anno)

# MINIMUM P VALUE IN DMR
min_DMR_p <- min_p_in_DMR(DMPs_anno, all_DMR_cpgs)


################################################################################
# MANHATTAN PLOT
################################################################################

# MANHATTAN PLOT
Manhattan_plot <- Manhattan(DMPs_anno, "FM_auto", signif_DMPs, all_DMR_cpgs)

ggsave("output/Figures/Manhattan_CaseCtrl_15CtrlPCs_2AncPCs.svg", Manhattan_plot, units = "cm", width = 13, height = 9, dpi = 300)
ggsave("output/Figures/Manhattan_CaseCtrl_15CtrlPCs_2AncPCs.pdf", Manhattan_plot, units = "cm", width = 13, height = 9, dpi = 300)
ggsave("output/Figures/Manhattan_CaseCtrl_15CtrlPCs_2AncPCs.png", Manhattan_plot, units = "cm", width = 13, height = 9, dpi = 400)


################################################################################
# ADJUST BETA VALUES
################################################################################

# SELECT M VALUES
# LOAD AND FILTER BETA VALUES
load(Mvalues_path)
Mvalues <- Mvalues[row.names(Mvalues) != "cg05575921", ] # remove cpg that is used to adjust for smoking
Mvalues_filt <- Mvalues[row.names(Mvalues) %in% DMPs_anno$cgID, colnames(Mvalues) %in% samples$Basename]


### PCAs ###
Ctrlprobe_PCAscores <- ctrl_probe_PCA(ctrl, samples)
cellcounts_PCAscores <- cell_type_PCA(samples, c("Epi", "Fib", "comb_ICs"))
snpsprobe_PCAscores <- ancestry_PCA(comb_SNPs, samples)

# adjust beta values
covar_matrix <- as.matrix(data.frame(Age = scale(as.numeric(samples$Age)), 
                                     cell_PC1 = cellcounts_PCAscores$PC1,
                                     cell_PC2 = cellcounts_PCAscores$PC2, 
                                     smoking = scale(as.numeric(samples$smoking)),
                                     Ctrl_PC1=Ctrlprobe_PCAscores$PC1, 
                                     Ctrl_PC2=Ctrlprobe_PCAscores$PC2, 
                                     Ctrl_PC3=Ctrlprobe_PCAscores$PC3, 
                                     Ctrl_PC4=Ctrlprobe_PCAscores$PC4, 
                                     Ctrl_PC5=Ctrlprobe_PCAscores$PC5,
                                     Ctrl_PC6=Ctrlprobe_PCAscores$PC6, 
                                     Ctrl_PC7=Ctrlprobe_PCAscores$PC7, 
                                     Ctrl_PC8=Ctrlprobe_PCAscores$PC8, 
                                     Ctrl_PC9=Ctrlprobe_PCAscores$PC9, 
                                     Ctrl_PC10=Ctrlprobe_PCAscores$PC10,
                                     Ctrl_PC11=Ctrlprobe_PCAscores$PC11,
                                     Ctrl_PC12=Ctrlprobe_PCAscores$PC12,
                                     Ctrl_PC13=Ctrlprobe_PCAscores$PC13,
                                     Ctrl_PC14=Ctrlprobe_PCAscores$PC14,
                                     Ctrl_PC15=Ctrlprobe_PCAscores$PC15,
                                     snps_PC1=snpsprobe_PCAscores$PC1, 
                                     snps_PC2=snpsprobe_PCAscores$PC2))

rm(Ctrlprobe_PCAscores, cellcounts_PCAscores, snpsprobe_PCAscores)

Mvalues_wobatch <- removeBatchEffect(Mvalues_filt, covariates = covar_matrix)
rm(covar_matrix)

beta_wobatch <- M2B(Mvalues_wobatch)
rm(Mvalues_wobatch)



################################################################################
# MANUAL ANNOTATION
################################################################################

dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED2")] <- "IDH3A"
dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED3")] <- "ARHGEF17"
dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED4")] <- "RGS12"
dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED5")] <- "BRSK2"
dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED6")] <- "MYO7A"
dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED7")] <- "TAPBP"
dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED8")] <- "UNANNOTATED2"
dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED9")] <- "UNANNOTATED3"
dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED10")] <- "LINC02580"
dmrff_sign$GeneString[which(dmrff_sign$GeneString == "UNANNOTATED11")] <- "PAX9"


################################################################################
# DMR VISUALISATION
################################################################################

dmrff_sign$avg_OCD <- NA
dmrff_sign$se_OCD <- NA
dmrff_sign$median_OCD <- NA
dmrff_sign$IQR_OCD <- NA

dmrff_sign$avg_CTRL <- NA
dmrff_sign$se_CTRL <- NA
dmrff_sign$median_CTRL <- NA
dmrff_sign$IQR_CTRL <- NA

dmrff_sign$top20th <- NA
dmrff_sign$bot20th <- NA
dmrff_sign$top20_OR <- NA
dmrff_sign$bottom20_OR <- NA
dmrff_sign$top20_OR_CI <- NA
dmrff_sign$bottom20_OR_CI <- NA
dmrff_sign$top20_OR_pvalue <- NA
dmrff_sign$bottom20_OR_pvalue <- NA



for (j in c(1:nrow(dmrff_sign))){
  
  chr <- dmrff_sign$chr[j]
  start <- dmrff_sign$start[j]
  end <- dmrff_sign$end[j]
  GeneString <- dmrff_sign$GeneString[j]
  GeneString2 <- gsub(";", "_", GeneString)
  DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
  beta_DMR <- beta_wobatch[rownames(beta_wobatch) %in% DMR_cpgs, ]
  
  
  # MAKE DF WITH CLINICAL INFO & AVERAGE DMR METHYLATION PER PATIENT
  DMR_categ <- add_avgDNAm(samples, beta_DMR, "Indiv_ID", c("Diagnosis"))
  
  ####################
  
  # CALCULATE ORs and RRs
  DMR_categ$Diagnosis <- as.vector(DMR_categ$Diagnosis)
  DMR_categ$Diagnosis <- factor(DMR_categ$Diagnosis, levels = c("CTRL", "OCD"))
  
  # CALCULATE QUANTILES
  top20_thresh <- quantile(DMR_categ$AvgDMRmethylation, 0.8)
  bottom20_thresh <- quantile(DMR_categ$AvgDMRmethylation, 0.2)
  dmrff_sign$top20th[j] <- top20_thresh
  dmrff_sign$bot20th[j] <- bottom20_thresh
  
  # HIGHLIGHT INDIVIDUALS IN TOP/BOTTOM QUANTILES
  DMR_categ$top20 <- DMR_categ$AvgDMRmethylation >= top20_thresh
  DMR_categ$bottom20 <- DMR_categ$AvgDMRmethylation <= bottom20_thresh
  
  # CALCULATE ORs
  table_top20 <- table(DMR_categ$Diagnosis, DMR_categ$top20)
  OR_top20 <- fisher.test(table_top20)
  
  table_bottom20 <- table(DMR_categ$Diagnosis, DMR_categ$bottom20)
  OR_bottom20 <- fisher.test(table_bottom20)
  
  dmrff_sign$top20_OR[j] <- OR_top20$estimate[1]
  dmrff_sign$bottom20_OR[j] <- OR_bottom20$estimate[1]
  dmrff_sign$top20_OR_CI[j] <- paste(round(OR_top20$conf.int[1:2],3), collapse = "-")
  dmrff_sign$bottom20_OR_CI[j] <- paste(round(OR_bottom20$conf.int[1:2],3), collapse = "-")
  dmrff_sign$top20_OR_pvalue[j] <- OR_top20$p.value[1]
  dmrff_sign$bottom20_OR_pvalue[j] <- OR_bottom20$p.value[1]
  
  
  if (OR_top20$p.value < 0.05){
    sign_top20 <- " *"
  } else {
    sign_top20 <- ""
  }
  
  if (OR_bottom20$p.value < 0.05){
    sign_bot20 <- " *"
  } else {
    sign_bot20 <- ""
  }
  
  # CALCULATE HIGHEST DENSITY VALUE
  dens_value <- density(DMR_categ$AvgDMRmethylation[which(DMR_categ$Diagnosis == "CTRL")])
  dens_value2 <- density(DMR_categ$AvgDMRmethylation[which(DMR_categ$Diagnosis == "OCD")])
  max_value <- max(c(dens_value$y, dens_value2$y))
  
  # PLOT  DENSITY PLOTS WITH ADDED ORs
  Dens_plot <- ggplot(DMR_categ, aes(x = AvgDMRmethylation, fill = Diagnosis)) + 
    geom_density(alpha = 0.5) + 
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    th + th_transparent +
    labs(x = "Average DNAm of CpGs in DMR", 
         y = "Density") +
    ggtitle(GeneString) +
    scale_y_continuous(expand=c(0,0), limits = c(0,max_value*1.3)) +
    scale_x_continuous(expand=c(0,0), limits = c(min(DMR_categ$AvgDMRmethylation),max(DMR_categ$AvgDMRmethylation))) +
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
    geom_segment(aes(x=top20_thresh, y = 0, xend=top20_thresh, yend = max_value), colour = "black") +
    annotate("text", x = top20_thresh, y=max_value, label = paste0("top 20%\nOR: ", round(OR_top20$estimate[1], 1), sign_top20), vjust = -0.25) +
    geom_segment(aes(x=bottom20_thresh, y = 0, xend=bottom20_thresh, yend = max_value), colour = "black") +
    annotate("text", x = bottom20_thresh, y=max_value, label = paste0("bottom 20%\nOR: ", round(OR_bottom20$estimate[1], 1), sign_bot20), vjust = -0.25)
  Dens_plot
  
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_DENS.svg"), Dens_plot, units = "cm", width = 13, height = 12, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_DENS.pdf"), Dens_plot, units = "cm", width = 13, height = 12, dpi = 400) 
  save(Dens_plot, max_value, OR_bottom20, OR_top20,bottom20_thresh, top20_thresh, DMR_categ,   
       file = paste0("output/Figures/DMRs/raw/DMR_CaseCtrl_", GeneString2, "_DENS.RData"))
  
  # CALCULATE THE MEANS AND STANDARD ERRORS FOR EACH TIME POINT
  df_summary <- as.data.frame(DMR_categ %>%
                                group_by(Diagnosis) %>%
                                summarize(Mean = mean(AvgDMRmethylation, na.rm = TRUE),
                                          se = sd(AvgDMRmethylation, na.rm = TRUE) / sqrt(n()),
                                          median = median(AvgDMRmethylation, na.rm = TRUE),
                                          IQR = IQR(AvgDMRmethylation, na.rm = TRUE)))
  rownames(df_summary) <- levels(DMR_categ$Diagnosis)
  
  # ADD TO OVERVIEW DF
  dmrff_sign$avg_OCD[j] <- round(df_summary["OCD","Mean"], 4)
  dmrff_sign$avg_CTRL[j] <- round(df_summary["CTRL","Mean"], 4)
  
  dmrff_sign$se_OCD[j] <- round(df_summary["OCD","se"], 4)
  dmrff_sign$se_CTRL[j] <- round(df_summary["CTRL","se"], 4)
  
  dmrff_sign$median_OCD[j] <- round(df_summary["OCD","median"], 4)
  dmrff_sign$median_CTRL[j] <- round(df_summary["CTRL","median"], 4)
  
  dmrff_sign$IQR_OCD[j] <- round(df_summary["OCD","IQR"], 4)
  dmrff_sign$IQR_CTRL[j] <- round(df_summary["CTRL","IQR"], 4)
  
  
  # CATEGORY PLOT
  Categ_plot <- ggplot(DMR_categ, aes(x = Diagnosis, y = AvgDMRmethylation, fill = Diagnosis)) +
    geom_violin(scale = "width", alpha = 0.3) + 
    geom_boxplot(width = 0.1, fill = "white") +
    geom_point(data = df_summary, aes(y = Mean, x = Diagnosis), color = "#14747b", size = 1) +
    th + th_transparent +
    ggtitle(GeneString) +
    labs(x = NULL, y = "Average DNAm of CpGs in DMR") + 
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a"))+
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_CATEG.pdf"), Categ_plot, units = "cm", width = 7, height = 10.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_CATEG.svg"), Categ_plot, units = "cm", width = 7, height = 10.5, dpi = 400)
  save(Categ_plot, df_summary, GeneString, file = paste0("output/Figures/DMRs/raw/DMR_CaseCtrl_", GeneString2, "_CATEG.RData"))
  
  
  # DMRFF PLOT
  dmrff_plot = dmrff.plot(region.chr = dmrff_sign$chr[j],
                          region.start = dmrff_sign$start[j], 
                          region.end = dmrff_sign$end[j], 
                          estimate = DMPs_anno$effectsize,
                          se = DMPs_anno$SE, 
                          chr = DMPs_anno$CHR,
                          pos = DMPs_anno$MAPINFO)
  
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_DMRFF.pdf"), dmrff_plot, units = "cm", width = 14, height = 10, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_DMRFF.svg"), dmrff_plot, units = "cm", width = 14, height = 10, dpi = 400)
  save(dmrff_plot, file = paste0("output/Figures/DMRs/raw/DMR_CaseCtrl_", GeneString2, "_DMRFF.RData"))
  
  
  # CPG PLOTS
  
  cpg_df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
  colnames(cpg_df) <- c("Basename", "Diagnosis", "cpg", "DNAm")
  
  DMR_cpgs_incl <- intersect(DMR_cpgs, rownames(beta_wobatch))
  
  for (cg in DMR_cpgs_incl){
    single_cg_df <- data.frame(Basename = samples$Basename, CaseCtrl = samples$Diagnosis, cpg = cg, DNAm = NA)
    
    for (i in c(1:nrow(single_cg_df))){
      Basename <- single_cg_df$Basename[i]
      betas_cpg <- beta_wobatch[cg, Basename]
      single_cg_df[i, "DNAm"] <- betas_cpg
    }
    
    cpg_df <- rbind(cpg_df, single_cg_df)
  }
  
  
  
  cpg_df$CaseCtrl <- factor(cpg_df$CaseCtrl, levels = c("CTRL", "OCD"))
  
  cpg_plot <- ggplot(cpg_df, aes(x=cpg, y = DNAm, fill = CaseCtrl)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    th + th_transparent + 
    ggtitle(GeneString) +
    labs(y = "DNA methylation") +
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(), 
          legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  cpg_plot
  
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_CPG.svg"), cpg_plot, units = "cm", width = 17, height = 13.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_CPG.pdf"), cpg_plot, units = "cm", width = 17, height = 13.5, dpi = 400)       
  save(cpg_plot, file = paste0("output/Figures/DMRs/raw/DMR_CaseCtrl_", GeneString2, "_CPG.RData"))
  
  
  median_data <- as.data.frame(cpg_df %>%
                                 group_by(cpg, CaseCtrl) %>%
                                 summarize(median = median(DNAm), .groups = "drop"))
  
  merged_cpg <- merge(x=median_data, y=anno, by.x = "cpg", by.y = "Name", all.x = TRUE)
  
  x_label <- paste0("position on ", chr)
  
  cpg_line_plot <- ggplot(merged_cpg, aes(x=MAPINFO, y = median, color = CaseCtrl, group = CaseCtrl)) +
    geom_point(size = 1) +
    geom_line() +
    labs(
      x = x_label,
      y = "median DNA methylation",
      color = NULL
    ) + th + th_transparent +
    ggtitle(GeneString) +
    scale_color_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    theme(legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  cpg_line_plot
  
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_TRACK.svg"), cpg_line_plot, units = "cm", width = 20, height = 13.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_TRACK.pdf"), cpg_line_plot, units = "cm", width = 20, height = 13.5, dpi = 400)       
  save(cpg_line_plot, GeneString, file = paste0("output/Figures/DMRs/raw/DMR_CaseCtrl_", GeneString2, "_TRACK.RData"))
  
  median_diff <- unique(merged_cpg[, c("cpg", "CHR", "MAPINFO")])
  median_diff$median_OCD_minus_CTRL <- NA
  
  for (i in c(1:nrow(median_diff))){
    cpg <- median_diff$cpg[i]
    OCD_value <- merged_cpg$median[which(merged_cpg$cpg == cpg & merged_cpg$CaseCtrl == "OCD")]
    CTRL_value <- merged_cpg$median[which(merged_cpg$cpg == cpg & merged_cpg$CaseCtrl == "CTRL")]
    median_diff$median_OCD_minus_CTRL[i] <- OCD_value - CTRL_value
  }
  
  cpg_line_plot_diff <- ggplot(median_diff, aes(x=MAPINFO, y = median_OCD_minus_CTRL)) +
    geom_point(size = 1, color = "#891a1a") +
    geom_line(color = "#891a1a") +
    labs(
      x = x_label,
      y = "difference in median DNAm (OCD-CTRL)",
      color = NULL
    ) + th + th_transparent +
    ggtitle(GeneString) +
    theme(legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    geom_hline(color = "#084081", yintercept = 0)
  cpg_line_plot_diff
  
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_TRACK_DIFF.svg"), cpg_line_plot_diff, units = "cm", width = 20, height = 10, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/DMR_CaseCtrl_", GeneString2, "_TRACK_DIFF.pdf"), cpg_line_plot_diff, units = "cm", width = 20, height = 10, dpi = 400)       
  save(cpg_line_plot_diff, GeneString, file = paste0("output/Figures/DMRs/raw/DMR_CaseCtrl_", GeneString2, "_TRACK_DIFF.RData"))
  
  
  
  
  
  
  
  
  Categ_plot <- Categ_plot + ggtitle("A")+ theme(plot.title = element_text(hjust = 0))
  Dens_plot <- Dens_plot + ggtitle("B")+ theme(plot.title = element_text(hjust = 0))
  #class(Dens_plot) <- class(Dens_plot)[class(Dens_plot) != "gg"]
  cpg_line_plot_diff <- cpg_line_plot_diff + ggtitle("C")+ theme(plot.title = element_text(hjust = 0))
  
  empty <- ggplot() + theme_void() #+ th + th_transparent
  
  design <- "
  ABBB
  CCCC
  DDDD
  "
  
  layout <- (Categ_plot + Dens_plot + empty + cpg_line_plot_diff) +
    plot_layout(design = design, heights = c(1, 0.1, 0.7)) +
    plot_annotation(title = GeneString,
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
  
  
  ggsave(paste0("output/Figures/DMRs/combined/DMR_CaseCtrl_", GeneString2, "_COMB.svg"), layout, units = "cm", width = 22, height = 24, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/combined/DMR_CaseCtrl_", GeneString2, "_COMB.pdf"), layout, units = "cm", width = 22, height = 24, dpi = 400)       
  
  rm(empty, layout)
  
}

dmrff_sign <- add_location(dmrff_sign, annotation)
dmrff_sign$avg_diff_OCD_Ctrl <- dmrff_sign$avg_OCD - dmrff_sign$avg_CTRL
dmrff_sign$median_diff_OCD_Ctrl <- dmrff_sign$median_OCD - dmrff_sign$median_CTRL

# EXPORT DMR TABLE WITH ADDED ORs
write_xlsx(dmrff_sign, path = "output/Tables/DMRs_CaseCtrl_15CtrlPCs_2AncPCs_extended.xlsx")
save(dmrff_sign, file = "output/RData/DMRs_CaseCtrl_15CtrlPCs_2AncPCs_extended.RData")






################################################################################
# ONLY UNMEDICATED
################################################################################

samples_unmed <- samples[which(samples$Psychoactive_medicine == "N" | (samples$Diagnosis == "CTRL" & is.na(samples$Psychoactive_medicine))), ]
table(samples_unmed$Diagnosis)
#CTRL  OCD 
#384  260 

beta_wobatch_unmed <- beta_wobatch[, colnames(beta_wobatch) %in% samples_unmed$Basename]
beta_unmed <- beta[, colnames(beta) %in% samples_unmed$Basename]

# prepare model
Ctrlprobe_PCAscores <- ctrl_probe_PCA(ctrl, samples_unmed)
cellcounts_PCAscores <- cell_type_PCA(samples_unmed, c("Epi", "Fib", "comb_ICs"))
snpsprobe_PCAscores <- ancestry_PCA(comb_SNPs, samples_unmed)

Sex = as.factor(samples_unmed$Sex)
Age = scale(as.numeric(samples_unmed$Age))
cell_PC1 = cellcounts_PCAscores$PC1
cell_PC2 = cellcounts_PCAscores$PC2 
smoking = scale(as.numeric(samples_unmed$smoking))
Ctrl_PC1=Ctrlprobe_PCAscores$PC1
Ctrl_PC2=Ctrlprobe_PCAscores$PC2
Ctrl_PC3=Ctrlprobe_PCAscores$PC3 
Ctrl_PC4=Ctrlprobe_PCAscores$PC4 
Ctrl_PC5=Ctrlprobe_PCAscores$PC5
Ctrl_PC6=Ctrlprobe_PCAscores$PC6 
Ctrl_PC7=Ctrlprobe_PCAscores$PC7 
Ctrl_PC8=Ctrlprobe_PCAscores$PC8 
Ctrl_PC9=Ctrlprobe_PCAscores$PC9 
Ctrl_PC10=Ctrlprobe_PCAscores$PC10
Ctrl_PC11=Ctrlprobe_PCAscores$PC11
Ctrl_PC12=Ctrlprobe_PCAscores$PC12
Ctrl_PC13=Ctrlprobe_PCAscores$PC13
Ctrl_PC14=Ctrlprobe_PCAscores$PC14
Ctrl_PC15=Ctrlprobe_PCAscores$PC15
snps_PC1=snpsprobe_PCAscores$PC1
snps_PC2=snpsprobe_PCAscores$PC2

unmedicated_result <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))
colnames(unmedicated_result) <- c("DMR", "Estimate", "SE", "t_value", "p_value", "median_diff_OCD_Ctrl", "all_median_diff")

for (j in c(1:nrow(dmrff_sign))){
  
  chr <- dmrff_sign$chr[j]
  start <- dmrff_sign$start[j]
  end <- dmrff_sign$end[j]
  GeneString <- dmrff_sign$GeneString[j]
  GeneString2 <- gsub(";", "_", GeneString)
  DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
  beta_DMR <- beta_unmed[rownames(beta_unmed) %in% DMR_cpgs, ]
  med_diff_old <- dmrff_sign$median_diff_OCD_Ctrl[j]
  
  
  # ASSOCIATION BETWEEN AVG DNAm AND CASE CONTROL STATUS
  DMR_categ <- add_avgDNAm(samples_unmed, beta_DMR, "Indiv_ID", c("Diagnosis"))
  DMR_categ$Avg_Mvalues <-  B2M(DMR_categ[, "AvgDMRmethylation"])
  
  avgDNAm <- DMR_categ$Avg_Mvalues
  CaseCtrl <- DMR_categ$Diagnosis
  
  lm_result <- summary(lm(avgDNAm ~ CaseCtrl + Sex + Age + cell_PC1 + cell_PC2 + smoking + 
       Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +
       Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
       Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
       snps_PC1 + snps_PC2
       ))
  
  result <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
  colnames(result) <- c("DMR", "Estimate", "SE", "t_value", "p_value", "median_diff_OCD_Ctrl", "all_median_diff")
  
  result[1, 1] <- GeneString2
  result[1, 2] <- lm_result$coefficients["CaseCtrlOCD", "Estimate"]
  result[1, 3] <- lm_result$coefficients["CaseCtrlOCD", "Std. Error"]
  result[1, 4] <- lm_result$coefficients["CaseCtrlOCD", "t value"]
  result[1, 5] <- lm_result$coefficients["CaseCtrlOCD", "Pr(>|t|)"]
  
  
  # MEDIAN DIFF OCD - CTRL
  beta_DMR2 <- beta_wobatch_unmed[rownames(beta_wobatch_unmed) %in% DMR_cpgs, ]
  DMR_categ2 <- add_avgDNAm(samples_unmed, beta_DMR2, "Indiv_ID", c("Diagnosis"))
  
  median_OCD <- median(DMR_categ2$AvgDMRmethylation[which(DMR_categ2$Diagnosis == "OCD")])
  median_CTRL <- median(DMR_categ2$AvgDMRmethylation[which(DMR_categ2$Diagnosis == "CTRL")])
  
  median_diff_OCD_Ctrl <- median_OCD - median_CTRL
  result[1, 6] <- median_diff_OCD_Ctrl
  
  result[1, 7] <- med_diff_old
  
  # COMBINE TABLES
  unmedicated_result <- rbind(unmedicated_result, result)
  
}

save(unmedicated_result, file = "output/RData/Sensitivity_unmedicated.RData")
write_xlsx(unmedicated_result, path = "output/Tables/Sensitivity_unmedicated.xlsx")


for (j in c(1:nrow(dmrff_sign))){
  
  chr <- dmrff_sign$chr[j]
  start <- dmrff_sign$start[j]
  end <- dmrff_sign$end[j]
  GeneString <- dmrff_sign$GeneString[j]
  GeneString2 <- gsub(";", "_", GeneString)
  DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
  beta_DMR <- beta_wobatch_unmed[rownames(beta_wobatch_unmed) %in% DMR_cpgs, ]
  
  
  # MAKE DF WITH CLINICAL INFO & AVERAGE DMR METHYLATION PER PATIENT
  DMR_categ <- add_avgDNAm(samples_unmed, beta_DMR, "Indiv_ID", c("Diagnosis"))
  
  ####################
  
  # CALCULATE ORs and RRs
  DMR_categ$Diagnosis <- as.vector(DMR_categ$Diagnosis)
  DMR_categ$Diagnosis <- factor(DMR_categ$Diagnosis, levels = c("CTRL", "OCD"))
  
  # CALCULATE QUANTILES
  top20_thresh <- quantile(DMR_categ$AvgDMRmethylation, 0.8)
  bottom20_thresh <- quantile(DMR_categ$AvgDMRmethylation, 0.2)
  
  # HIGHLIGHT INDIVIDUALS IN TOP/BOTTOM QUANTILES
  DMR_categ$top20 <- DMR_categ$AvgDMRmethylation >= top20_thresh
  DMR_categ$bottom20 <- DMR_categ$AvgDMRmethylation <= bottom20_thresh
  
  # CALCULATE ORs
  table_top20 <- table(DMR_categ$Diagnosis, DMR_categ$top20)
  OR_top20 <- fisher.test(table_top20)
  
  table_bottom20 <- table(DMR_categ$Diagnosis, DMR_categ$bottom20)
  OR_bottom20 <- fisher.test(table_bottom20)
  
  if (OR_top20$p.value < 0.05){
    sign_top20 <- " *"
  } else {
    sign_top20 <- ""
  }
  
  if (OR_bottom20$p.value < 0.05){
    sign_bot20 <- " *"
  } else {
    sign_bot20 <- ""
  }
  
  # CALCULATE HIGHEST DENSITY VALUE
  dens_value <- density(DMR_categ$AvgDMRmethylation[which(DMR_categ$Diagnosis == "CTRL")])
  dens_value2 <- density(DMR_categ$AvgDMRmethylation[which(DMR_categ$Diagnosis == "OCD")])
  max_value <- max(c(dens_value$y, dens_value2$y))
  
  # PLOT  DENSITY PLOTS WITH ADDED ORs
  Dens_plot <- ggplot(DMR_categ, aes(x = AvgDMRmethylation, fill = Diagnosis)) + 
    geom_density(alpha = 0.5) + 
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    th + th_transparent +
    labs(x = "Average DNAm of CpGs in DMR", 
         y = "Density") +
    ggtitle(GeneString) +
    scale_y_continuous(expand=c(0,0), limits = c(0,max_value*1.3)) +
    scale_x_continuous(expand=c(0,0), limits = c(min(DMR_categ$AvgDMRmethylation),max(DMR_categ$AvgDMRmethylation))) +
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
    geom_segment(aes(x=top20_thresh, y = 0, xend=top20_thresh, yend = max_value), colour = "black") +
    annotate("text", x = top20_thresh, y=max_value, label = paste0("top 20%\nOR: ", round(OR_top20$estimate[1], 1), sign_top20), vjust = -0.25) +
    geom_segment(aes(x=bottom20_thresh, y = 0, xend=bottom20_thresh, yend = max_value), colour = "black") +
    annotate("text", x = bottom20_thresh, y=max_value, label = paste0("bottom 20%\nOR: ", round(OR_bottom20$estimate[1], 1), sign_bot20), vjust = -0.25)
  Dens_plot
  
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_DENS_UNMEDICATED.svg"), Dens_plot, units = "cm", width = 13, height = 12, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_DENS_UNMEDICATED.pdf"), Dens_plot, units = "cm", width = 13, height = 12, dpi = 400)       
  
  # CALCULATE THE MEANS AND STANDARD ERRORS FOR EACH TIME POINT
  df_summary <- as.data.frame(DMR_categ %>%
                                group_by(Diagnosis) %>%
                                summarize(Mean = mean(AvgDMRmethylation, na.rm = TRUE),
                                          se = sd(AvgDMRmethylation, na.rm = TRUE) / sqrt(n()),
                                          median = median(AvgDMRmethylation, na.rm = TRUE),
                                          IQR = IQR(AvgDMRmethylation, na.rm = TRUE)))
  rownames(df_summary) <- levels(DMR_categ$Diagnosis)
  
  # CATEGORY PLOT
  Categ_plot <- ggplot(DMR_categ, aes(x = Diagnosis, y = AvgDMRmethylation, fill = Diagnosis)) +
    geom_violin(scale = "width", alpha = 0.3) + 
    geom_boxplot(width = 0.1, fill = "white") +
    geom_point(data = df_summary, aes(y = Mean, x = Diagnosis), color = "#14747b", size = 1) +
    th + th_transparent +
    ggtitle(GeneString) +
    labs(x = NULL, y = "Average DNAm of CpGs in DMR") + 
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a"))+
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
  Categ_plot
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_CATEG.pdf"), Categ_plot, units = "cm", width = 7, height = 10.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_CATEG.svg"), Categ_plot, units = "cm", width = 7, height = 10.5, dpi = 400)
  
  # CPG PLOTS
  
  cpg_df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
  colnames(cpg_df) <- c("Basename", "Diagnosis", "cpg", "DNAm")
  
  DMR_cpgs_incl <- intersect(DMR_cpgs, rownames(beta_wobatch_unmed))
  
  for (cg in DMR_cpgs_incl){
    single_cg_df <- data.frame(Basename = samples_unmed$Basename, CaseCtrl = samples_unmed$Diagnosis, cpg = cg, DNAm = NA)
    
    for (i in c(1:nrow(single_cg_df))){
      Basename <- single_cg_df$Basename[i]
      betas_cpg <- beta_wobatch_unmed[cg, Basename]
      single_cg_df[i, "DNAm"] <- betas_cpg
    }
    
    cpg_df <- rbind(cpg_df, single_cg_df)
  }
  
  
  
  cpg_df$CaseCtrl <- factor(cpg_df$CaseCtrl, levels = c("CTRL", "OCD"))
  
  cpg_plot <- ggplot(cpg_df, aes(x=cpg, y = DNAm, fill = CaseCtrl)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    th + th_transparent + 
    ggtitle(GeneString) +
    labs(y = "DNA methylation") +
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(), 
          legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  cpg_plot
  
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_CPG.svg"), cpg_plot, units = "cm", width = 17, height = 13.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_CPG.pdf"), cpg_plot, units = "cm", width = 17, height = 13.5, dpi = 400)       
  
  
  
  median_data <- as.data.frame(cpg_df %>%
                                 group_by(cpg, CaseCtrl) %>%
                                 summarize(median = median(DNAm), .groups = "drop"))
  
  merged_cpg <- merge(x=median_data, y=anno, by.x = "cpg", by.y = "Name", all.x = TRUE)
  
  x_label <- paste0("position on ", chr)
  
  cpg_line_plot <- ggplot(merged_cpg, aes(x=MAPINFO, y = median, color = CaseCtrl, group = CaseCtrl)) +
    geom_point(size = 1) +
    geom_line() +
    labs(
      x = x_label,
      y = "median DNA methylation",
      color = NULL
    ) + th + th_transparent +
    ggtitle(GeneString) +
    scale_color_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    theme(legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  cpg_line_plot
  
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_TRACK.svg"), cpg_line_plot, units = "cm", width = 20, height = 13.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_TRACK.pdf"), cpg_line_plot, units = "cm", width = 20, height = 13.5, dpi = 400)       
  
  
  median_diff <- unique(merged_cpg[, c("cpg", "CHR", "MAPINFO")])
  median_diff$median_OCD_minus_CTRL <- NA
  
  for (i in c(1:nrow(median_diff))){
    cpg <- median_diff$cpg[i]
    OCD_value <- merged_cpg$median[which(merged_cpg$cpg == cpg & merged_cpg$CaseCtrl == "OCD")]
    CTRL_value <- merged_cpg$median[which(merged_cpg$cpg == cpg & merged_cpg$CaseCtrl == "CTRL")]
    median_diff$median_OCD_minus_CTRL[i] <- OCD_value - CTRL_value
  }
  
  cpg_line_plot_diff <- ggplot(median_diff, aes(x=MAPINFO, y = median_OCD_minus_CTRL)) +
    geom_point(size = 1, color = "#891a1a") +
    geom_line(color = "#891a1a") +
    labs(
      x = x_label,
      y = "difference in median DNAm (OCD-CTRL)",
      color = NULL
    ) + th + th_transparent +
    ggtitle(GeneString) +
    theme(legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    geom_hline(color = "#084081", yintercept = 0)
  cpg_line_plot_diff
  
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_TRACK_DIFF.svg"), cpg_line_plot_diff, units = "cm", width = 20, height = 10, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/unmedicated/DMR_CaseCtrl_", GeneString2, "_TRACK_DIFF.pdf"), cpg_line_plot_diff, units = "cm", width = 20, height = 10, dpi = 400)       
  
  
  
  Categ_plot <- Categ_plot + ggtitle("A")+ theme(plot.title = element_text(hjust = 0))
  Dens_plot <- Dens_plot + ggtitle("B")+ theme(plot.title = element_text(hjust = 0))
  #class(Dens_plot) <- class(Dens_plot)[class(Dens_plot) != "gg"]
  cpg_line_plot_diff <- cpg_line_plot_diff + ggtitle("C")+ theme(plot.title = element_text(hjust = 0))
  
  empty <- ggplot() + theme_void() #+ th + th_transparent
  
  design <- "
  ABBB
  CCCC
  DDDD
  "
  
  layout <- (Categ_plot + Dens_plot + empty + cpg_line_plot_diff) +
    plot_layout(design = design, heights = c(1, 0.1, 0.7)) +
    plot_annotation(title = GeneString,
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
  
  
  ggsave(paste0("output/Figures/DMRs/unmedicated/combined/DMR_CaseCtrl_", GeneString2, "_COMB.svg"), layout, units = "cm", width = 22, height = 24, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/unmedicated/combined/DMR_CaseCtrl_", GeneString2, "_COMB.pdf"), layout, units = "cm", width = 22, height = 24, dpi = 400)       
  
  rm(empty, layout)
  
}




################################################################################
# ONLY EUROPEAN
################################################################################

samples_eur <- samples[which(samples$Ancestry == "EUR"), ]
table(samples_eur$Diagnosis)
#CTRL  OCD 
#132  370 


beta_wobatch_eur <- beta_wobatch[, colnames(beta_wobatch) %in% samples_eur$Basename]
beta_eur <- beta[, colnames(beta) %in% samples_eur$Basename]

# prepare model
Ctrlprobe_PCAscores <- ctrl_probe_PCA(ctrl, samples_eur)
cellcounts_PCAscores <- cell_type_PCA(samples_eur, c("Epi", "Fib", "comb_ICs"))
snpsprobe_PCAscores <- ancestry_PCA(comb_SNPs, samples_eur)

Sex = as.factor(samples_eur$Sex)
Age = scale(as.numeric(samples_eur$Age))
cell_PC1 = cellcounts_PCAscores$PC1
cell_PC2 = cellcounts_PCAscores$PC2 
smoking = scale(as.numeric(samples_eur$smoking))
Ctrl_PC1=Ctrlprobe_PCAscores$PC1
Ctrl_PC2=Ctrlprobe_PCAscores$PC2
Ctrl_PC3=Ctrlprobe_PCAscores$PC3 
Ctrl_PC4=Ctrlprobe_PCAscores$PC4 
Ctrl_PC5=Ctrlprobe_PCAscores$PC5
Ctrl_PC6=Ctrlprobe_PCAscores$PC6 
Ctrl_PC7=Ctrlprobe_PCAscores$PC7 
Ctrl_PC8=Ctrlprobe_PCAscores$PC8 
Ctrl_PC9=Ctrlprobe_PCAscores$PC9 
Ctrl_PC10=Ctrlprobe_PCAscores$PC10
Ctrl_PC11=Ctrlprobe_PCAscores$PC11
Ctrl_PC12=Ctrlprobe_PCAscores$PC12
Ctrl_PC13=Ctrlprobe_PCAscores$PC13
Ctrl_PC14=Ctrlprobe_PCAscores$PC14
Ctrl_PC15=Ctrlprobe_PCAscores$PC15
snps_PC1=snpsprobe_PCAscores$PC1
snps_PC2=snpsprobe_PCAscores$PC2

eur_result <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))
colnames(eur_result) <- c("DMR", "Estimate", "SE", "t_value", "p_value", "median_diff_OCD_Ctrl", "all_median_diff")

for (j in c(1:nrow(dmrff_sign))){
  
  chr <- dmrff_sign$chr[j]
  start <- dmrff_sign$start[j]
  end <- dmrff_sign$end[j]
  GeneString <- dmrff_sign$GeneString[j]
  GeneString2 <- gsub(";", "_", GeneString)
  DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
  beta_DMR <- beta_eur[rownames(beta_eur) %in% DMR_cpgs, ]
  med_diff_old <- dmrff_sign$median_diff_OCD_Ctrl[j]
  

  # ASSOCIATION BETWEEN AVG DNAm AND CASE CONTROL STATUS
  DMR_categ <- add_avgDNAm(samples_eur, beta_DMR, "Indiv_ID", c("Diagnosis"))
  DMR_categ$Avg_Mvalues <-  B2M(DMR_categ[, "AvgDMRmethylation"])
  
  avgDNAm <- DMR_categ$Avg_Mvalues
  CaseCtrl <- DMR_categ$Diagnosis
  
  lm_result <- summary(lm(avgDNAm ~ CaseCtrl + Sex + Age + cell_PC1 + cell_PC2 + smoking + 
                            Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +
                            Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                            Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                            snps_PC1 + snps_PC2
  ))
  
  result <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
  colnames(result) <- c("DMR", "Estimate", "SE", "t_value", "p_value", "median_diff_OCD_Ctrl", "all_median_diff")
  
  result[1, 1] <- GeneString2
  result[1, 2] <- lm_result$coefficients["CaseCtrlOCD", "Estimate"]
  result[1, 3] <- lm_result$coefficients["CaseCtrlOCD", "Std. Error"]
  result[1, 4] <- lm_result$coefficients["CaseCtrlOCD", "t value"]
  result[1, 5] <- lm_result$coefficients["CaseCtrlOCD", "Pr(>|t|)"]
  
  
  # MEDIAN DIFF OCD - CTRL
  beta_DMR2 <- beta_wobatch_eur[rownames(beta_wobatch_eur) %in% DMR_cpgs, ]
  DMR_categ2 <- add_avgDNAm(samples_eur, beta_DMR2, "Indiv_ID", c("Diagnosis"))
  
  median_OCD <- median(DMR_categ2$AvgDMRmethylation[which(DMR_categ2$Diagnosis == "OCD")])
  median_CTRL <- median(DMR_categ2$AvgDMRmethylation[which(DMR_categ2$Diagnosis == "CTRL")])
  
  median_diff_OCD_Ctrl <- median_OCD - median_CTRL
  result[1, 6] <- median_diff_OCD_Ctrl
  
  result[1, 7] <- med_diff_old
  
  # COMBINE TABLES
  eur_result <- rbind(eur_result, result)
  
}

save(eur_result, file = "output/RData/Sensitivity_european.RData")
write_xlsx(eur_result, path = "output/Tables/Sensitivity_european.xlsx")



for (j in c(1:nrow(dmrff_sign))){
  
  chr <- dmrff_sign$chr[j]
  start <- dmrff_sign$start[j]
  end <- dmrff_sign$end[j]
  GeneString <- dmrff_sign$GeneString[j]
  GeneString2 <- gsub(";", "_", GeneString)
  DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
  beta_DMR <- beta_wobatch_eur[rownames(beta_wobatch_eur) %in% DMR_cpgs, ]
  
  
  # MAKE DF WITH CLINICAL INFO & AVERAGE DMR METHYLATION PER PATIENT
  DMR_categ <- add_avgDNAm(samples_eur, beta_DMR, "Indiv_ID", c("Diagnosis"))
  
  ####################
  
  # CALCULATE ORs and RRs
  DMR_categ$Diagnosis <- as.vector(DMR_categ$Diagnosis)
  DMR_categ$Diagnosis <- factor(DMR_categ$Diagnosis, levels = c("CTRL", "OCD"))
  
  # CALCULATE QUANTILES
  top20_thresh <- quantile(DMR_categ$AvgDMRmethylation, 0.8)
  bottom20_thresh <- quantile(DMR_categ$AvgDMRmethylation, 0.2)
  
  # HIGHLIGHT INDIVIDUALS IN TOP/BOTTOM QUANTILES
  DMR_categ$top20 <- DMR_categ$AvgDMRmethylation >= top20_thresh
  DMR_categ$bottom20 <- DMR_categ$AvgDMRmethylation <= bottom20_thresh
  
  # CALCULATE ORs
  table_top20 <- table(DMR_categ$Diagnosis, DMR_categ$top20)
  OR_top20 <- fisher.test(table_top20)
  
  table_bottom20 <- table(DMR_categ$Diagnosis, DMR_categ$bottom20)
  OR_bottom20 <- fisher.test(table_bottom20)
  
  if (OR_top20$p.value < 0.05){
    sign_top20 <- " *"
  } else {
    sign_top20 <- ""
  }
  
  if (OR_bottom20$p.value < 0.05){
    sign_bot20 <- " *"
  } else {
    sign_bot20 <- ""
  }
  
  # CALCULATE HIGHEST DENSITY VALUE
  dens_value <- density(DMR_categ$AvgDMRmethylation[which(DMR_categ$Diagnosis == "CTRL")])
  dens_value2 <- density(DMR_categ$AvgDMRmethylation[which(DMR_categ$Diagnosis == "OCD")])
  max_value <- max(c(dens_value$y, dens_value2$y))
  
  # PLOT  DENSITY PLOTS WITH ADDED ORs
  Dens_plot <- ggplot(DMR_categ, aes(x = AvgDMRmethylation, fill = Diagnosis)) + 
    geom_density(alpha = 0.5) + 
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    th + th_transparent +
    labs(x = "Average DNAm of CpGs in DMR", 
         y = "Density") +
    ggtitle(GeneString) +
    scale_y_continuous(expand=c(0,0), limits = c(0,max_value*1.3)) +
    scale_x_continuous(expand=c(0,0), limits = c(min(DMR_categ$AvgDMRmethylation),max(DMR_categ$AvgDMRmethylation))) +
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
    geom_segment(aes(x=top20_thresh, y = 0, xend=top20_thresh, yend = max_value), colour = "black") +
    annotate("text", x = top20_thresh, y=max_value, label = paste0("top 20%\nOR: ", round(OR_top20$estimate[1], 1), sign_top20), vjust = -0.25) +
    geom_segment(aes(x=bottom20_thresh, y = 0, xend=bottom20_thresh, yend = max_value), colour = "black") +
    annotate("text", x = bottom20_thresh, y=max_value, label = paste0("bottom 20%\nOR: ", round(OR_bottom20$estimate[1], 1), sign_bot20), vjust = -0.25)
  Dens_plot
  
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_DENS_european.svg"), Dens_plot, units = "cm", width = 13, height = 12, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_DENS_european.pdf"), Dens_plot, units = "cm", width = 13, height = 12, dpi = 400)       
  
  # CALCULATE THE MEANS AND STANDARD ERRORS FOR EACH TIME POINT
  df_summary <- as.data.frame(DMR_categ %>%
                                group_by(Diagnosis) %>%
                                summarize(Mean = mean(AvgDMRmethylation, na.rm = TRUE),
                                          se = sd(AvgDMRmethylation, na.rm = TRUE) / sqrt(n()),
                                          median = median(AvgDMRmethylation, na.rm = TRUE),
                                          IQR = IQR(AvgDMRmethylation, na.rm = TRUE)))
  rownames(df_summary) <- levels(DMR_categ$Diagnosis)
  
  # CATEGORY PLOT
  Categ_plot <- ggplot(DMR_categ, aes(x = Diagnosis, y = AvgDMRmethylation, fill = Diagnosis)) +
    geom_violin(scale = "width", alpha = 0.3) + 
    geom_boxplot(width = 0.1, fill = "white") +
    geom_point(data = df_summary, aes(y = Mean, x = Diagnosis), color = "#14747b", size = 1) +
    th + th_transparent +
    ggtitle(GeneString) +
    labs(x = NULL, y = "Average DNAm of CpGs in DMR") + 
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a"))+
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
  Categ_plot
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_CATEG.pdf"), Categ_plot, units = "cm", width = 7, height = 10.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_CATEG.svg"), Categ_plot, units = "cm", width = 7, height = 10.5, dpi = 400)
  
  # CPG PLOTS
  
  cpg_df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
  colnames(cpg_df) <- c("Basename", "Diagnosis", "cpg", "DNAm")
  
  DMR_cpgs_incl <- intersect(DMR_cpgs, rownames(beta_wobatch_unmed))
  
  for (cg in DMR_cpgs_incl){
    single_cg_df <- data.frame(Basename = samples_eur$Basename, CaseCtrl = samples_eur$Diagnosis, cpg = cg, DNAm = NA)
    
    for (i in c(1:nrow(single_cg_df))){
      Basename <- single_cg_df$Basename[i]
      betas_cpg <- beta_wobatch_unmed[cg, Basename]
      single_cg_df[i, "DNAm"] <- betas_cpg
    }
    
    cpg_df <- rbind(cpg_df, single_cg_df)
  }
  
  
  
  cpg_df$CaseCtrl <- factor(cpg_df$CaseCtrl, levels = c("CTRL", "OCD"))
  
  cpg_plot <- ggplot(cpg_df, aes(x=cpg, y = DNAm, fill = CaseCtrl)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    th + th_transparent + 
    ggtitle(GeneString) +
    labs(y = "DNA methylation") +
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(), 
          legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  cpg_plot
  
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_CPG.svg"), cpg_plot, units = "cm", width = 17, height = 13.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_CPG.pdf"), cpg_plot, units = "cm", width = 17, height = 13.5, dpi = 400)       
  
  
  
  median_data <- as.data.frame(cpg_df %>%
                                 group_by(cpg, CaseCtrl) %>%
                                 summarize(median = median(DNAm), .groups = "drop"))
  
  merged_cpg <- merge(x=median_data, y=anno, by.x = "cpg", by.y = "Name", all.x = TRUE)
  
  x_label <- paste0("position on ", chr)
  
  cpg_line_plot <- ggplot(merged_cpg, aes(x=MAPINFO, y = median, color = CaseCtrl, group = CaseCtrl)) +
    geom_point(size = 1) +
    geom_line() +
    labs(
      x = x_label,
      y = "median DNA methylation",
      color = NULL
    ) + th + th_transparent +
    ggtitle(GeneString) +
    scale_color_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    theme(legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  cpg_line_plot
  
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_TRACK.svg"), cpg_line_plot, units = "cm", width = 20, height = 13.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_TRACK.pdf"), cpg_line_plot, units = "cm", width = 20, height = 13.5, dpi = 400)       
  
  
  median_diff <- unique(merged_cpg[, c("cpg", "CHR", "MAPINFO")])
  median_diff$median_OCD_minus_CTRL <- NA
  
  for (i in c(1:nrow(median_diff))){
    cpg <- median_diff$cpg[i]
    OCD_value <- merged_cpg$median[which(merged_cpg$cpg == cpg & merged_cpg$CaseCtrl == "OCD")]
    CTRL_value <- merged_cpg$median[which(merged_cpg$cpg == cpg & merged_cpg$CaseCtrl == "CTRL")]
    median_diff$median_OCD_minus_CTRL[i] <- OCD_value - CTRL_value
  }
  
  cpg_line_plot_diff <- ggplot(median_diff, aes(x=MAPINFO, y = median_OCD_minus_CTRL)) +
    geom_point(size = 1, color = "#891a1a") +
    geom_line(color = "#891a1a") +
    labs(
      x = x_label,
      y = "difference in median DNAm (OCD-CTRL)",
      color = NULL
    ) + th + th_transparent +
    ggtitle(GeneString) +
    theme(legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    geom_hline(color = "#084081", yintercept = 0)
  cpg_line_plot_diff
  
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_TRACK_DIFF.svg"), cpg_line_plot_diff, units = "cm", width = 20, height = 10, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/european/DMR_CaseCtrl_", GeneString2, "_TRACK_DIFF.pdf"), cpg_line_plot_diff, units = "cm", width = 20, height = 10, dpi = 400)       
  
  
  
  
  
  
  
  
  
  Categ_plot <- Categ_plot + ggtitle("A")+ theme(plot.title = element_text(hjust = 0))
  Dens_plot <- Dens_plot + ggtitle("B")+ theme(plot.title = element_text(hjust = 0))
  #class(Dens_plot) <- class(Dens_plot)[class(Dens_plot) != "gg"]
  cpg_line_plot_diff <- cpg_line_plot_diff + ggtitle("C")+ theme(plot.title = element_text(hjust = 0))
  
  empty <- ggplot() + theme_void() #+ th + th_transparent
  
  design <- "
  ABBB
  CCCC
  DDDD
  "
  
  layout <- (Categ_plot + Dens_plot + empty + cpg_line_plot_diff) +
    plot_layout(design = design, heights = c(1, 0.1, 0.7)) +
    plot_annotation(title = GeneString,
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
  
  
  ggsave(paste0("output/Figures/DMRs/european/combined/DMR_CaseCtrl_", GeneString2, "_COMB.svg"), layout, units = "cm", width = 22, height = 24, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/european/combined/DMR_CaseCtrl_", GeneString2, "_COMB.pdf"), layout, units = "cm", width = 22, height = 24, dpi = 400)       
  
  rm(empty, layout)
  
}







################################################################################
# NO COMORBIDITY
################################################################################

samples_noComorb <- samples[which(samples$Diagnosis == "CTRL" | (samples$Diagnosis == "OCD" & samples$Comorbidity == "N")), ]
table(samples_noComorb$Diagnosis)
#CTRL  OCD 
#384  227 


beta_wobatch_noComorb <- beta_wobatch[, colnames(beta_wobatch) %in% samples_noComorb$Basename]
beta_noComorb <- beta[, colnames(beta) %in% samples_noComorb$Basename]

# prepare model
Ctrlprobe_PCAscores <- ctrl_probe_PCA(ctrl, samples_noComorb)
cellcounts_PCAscores <- cell_type_PCA(samples_noComorb, c("Epi", "Fib", "comb_ICs"))
snpsprobe_PCAscores <- ancestry_PCA(comb_SNPs, samples_noComorb)

Sex = as.factor(samples_noComorb$Sex)
Age = scale(as.numeric(samples_noComorb$Age))
cell_PC1 = cellcounts_PCAscores$PC1
cell_PC2 = cellcounts_PCAscores$PC2 
smoking = scale(as.numeric(samples_noComorb$smoking))
Ctrl_PC1=Ctrlprobe_PCAscores$PC1
Ctrl_PC2=Ctrlprobe_PCAscores$PC2
Ctrl_PC3=Ctrlprobe_PCAscores$PC3 
Ctrl_PC4=Ctrlprobe_PCAscores$PC4 
Ctrl_PC5=Ctrlprobe_PCAscores$PC5
Ctrl_PC6=Ctrlprobe_PCAscores$PC6 
Ctrl_PC7=Ctrlprobe_PCAscores$PC7 
Ctrl_PC8=Ctrlprobe_PCAscores$PC8 
Ctrl_PC9=Ctrlprobe_PCAscores$PC9 
Ctrl_PC10=Ctrlprobe_PCAscores$PC10
Ctrl_PC11=Ctrlprobe_PCAscores$PC11
Ctrl_PC12=Ctrlprobe_PCAscores$PC12
Ctrl_PC13=Ctrlprobe_PCAscores$PC13
Ctrl_PC14=Ctrlprobe_PCAscores$PC14
Ctrl_PC15=Ctrlprobe_PCAscores$PC15
snps_PC1=snpsprobe_PCAscores$PC1
snps_PC2=snpsprobe_PCAscores$PC2

noComorb_result <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))
colnames(noComorb_result) <- c("DMR", "Estimate", "SE", "t_value", "p_value", "median_diff_OCD_Ctrl", "all_median_diff")

for (j in c(1:nrow(dmrff_sign))){
  
  chr <- dmrff_sign$chr[j]
  start <- dmrff_sign$start[j]
  end <- dmrff_sign$end[j]
  GeneString <- dmrff_sign$GeneString[j]
  GeneString2 <- gsub(";", "_", GeneString)
  DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
  beta_DMR <- beta_noComorb[rownames(beta_noComorb) %in% DMR_cpgs, ]
  med_diff_old <- dmrff_sign$median_diff_OCD_Ctrl[j]
  
  
  
  # ASSOCIATION BETWEEN AVG DNAm AND CASE CONTROL STATUS
  DMR_categ <- add_avgDNAm(samples_noComorb, beta_DMR, "Indiv_ID", c("Diagnosis"))
  DMR_categ$Avg_Mvalues <-  B2M(DMR_categ[, "AvgDMRmethylation"])
  
  avgDNAm <- DMR_categ$Avg_Mvalues
  CaseCtrl <- DMR_categ$Diagnosis
  
  lm_result <- summary(lm(avgDNAm ~ CaseCtrl + Sex + Age + cell_PC1 + cell_PC2 + smoking + 
                            Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC4 + Ctrl_PC5 +
                            Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                            Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                            snps_PC1 + snps_PC2
  ))
  
  result <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 7))
  colnames(result) <- c("DMR", "Estimate", "SE", "t_value", "p_value", "median_diff_OCD_Ctrl", "all_median_diff")
  
  result[1, 1] <- GeneString2
  result[1, 2] <- lm_result$coefficients["CaseCtrlOCD", "Estimate"]
  result[1, 3] <- lm_result$coefficients["CaseCtrlOCD", "Std. Error"]
  result[1, 4] <- lm_result$coefficients["CaseCtrlOCD", "t value"]
  result[1, 5] <- lm_result$coefficients["CaseCtrlOCD", "Pr(>|t|)"]
  
  
  # MEDIAN DIFF OCD - CTRL
  beta_DMR2 <- beta_wobatch_noComorb[rownames(beta_wobatch_noComorb) %in% DMR_cpgs, ]
  DMR_categ2 <- add_avgDNAm(samples_noComorb, beta_DMR2, "Indiv_ID", c("Diagnosis"))
  
  median_OCD <- median(DMR_categ2$AvgDMRmethylation[which(DMR_categ2$Diagnosis == "OCD")])
  median_CTRL <- median(DMR_categ2$AvgDMRmethylation[which(DMR_categ2$Diagnosis == "CTRL")])
  
  median_diff_OCD_Ctrl <- median_OCD - median_CTRL
  result[1, 6] <- median_diff_OCD_Ctrl
  
  result[1, 7] <- med_diff_old
  
  # COMBINE TABLES
  noComorb_result <- rbind(noComorb_result, result)
  
}

save(noComorb_result, file = "output/RData/Sensitivity_noComorb.RData")
write_xlsx(noComorb_result, path = "output/Tables/Sensitivity_noComorb.xlsx")

for (j in c(1:nrow(dmrff_sign))){
  
  chr <- dmrff_sign$chr[j]
  start <- dmrff_sign$start[j]
  end <- dmrff_sign$end[j]
  GeneString <- dmrff_sign$GeneString[j]
  GeneString2 <- gsub(";", "_", GeneString)
  DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
  beta_DMR <- beta_wobatch_noComorb[rownames(beta_wobatch_noComorb) %in% DMR_cpgs, ]
  
  
  # MAKE DF WITH CLINICAL INFO & AVERAGE DMR METHYLATION PER PATIENT
  DMR_categ <- add_avgDNAm(samples_noComorb, beta_DMR, "Indiv_ID", c("Diagnosis"))
  
  ####################
  
  # CALCULATE ORs and RRs
  DMR_categ$Diagnosis <- as.vector(DMR_categ$Diagnosis)
  DMR_categ$Diagnosis <- factor(DMR_categ$Diagnosis, levels = c("CTRL", "OCD"))
  
  # CALCULATE QUANTILES
  top20_thresh <- quantile(DMR_categ$AvgDMRmethylation, 0.8)
  bottom20_thresh <- quantile(DMR_categ$AvgDMRmethylation, 0.2)
  
  # HIGHLIGHT INDIVIDUALS IN TOP/BOTTOM QUANTILES
  DMR_categ$top20 <- DMR_categ$AvgDMRmethylation >= top20_thresh
  DMR_categ$bottom20 <- DMR_categ$AvgDMRmethylation <= bottom20_thresh
  
  # CALCULATE ORs
  table_top20 <- table(DMR_categ$Diagnosis, DMR_categ$top20)
  OR_top20 <- fisher.test(table_top20)
  
  table_bottom20 <- table(DMR_categ$Diagnosis, DMR_categ$bottom20)
  OR_bottom20 <- fisher.test(table_bottom20)
  
  if (OR_top20$p.value < 0.05){
    sign_top20 <- " *"
  } else {
    sign_top20 <- ""
  }
  
  if (OR_bottom20$p.value < 0.05){
    sign_bot20 <- " *"
  } else {
    sign_bot20 <- ""
  }
  
  # CALCULATE HIGHEST DENSITY VALUE
  dens_value <- density(DMR_categ$AvgDMRmethylation[which(DMR_categ$Diagnosis == "CTRL")])
  dens_value2 <- density(DMR_categ$AvgDMRmethylation[which(DMR_categ$Diagnosis == "OCD")])
  max_value <- max(c(dens_value$y, dens_value2$y))
  
  # PLOT  DENSITY PLOTS WITH ADDED ORs
  Dens_plot <- ggplot(DMR_categ, aes(x = AvgDMRmethylation, fill = Diagnosis)) + 
    geom_density(alpha = 0.5) + 
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    th + th_transparent +
    labs(x = "Average DNAm of CpGs in DMR", 
         y = "Density") +
    ggtitle(GeneString) +
    scale_y_continuous(expand=c(0,0), limits = c(0,max_value*1.3)) +
    scale_x_continuous(expand=c(0,0), limits = c(min(DMR_categ$AvgDMRmethylation),max(DMR_categ$AvgDMRmethylation))) +
    theme(legend.position = "top", legend.title = element_blank(), legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5)) +
    geom_segment(aes(x=top20_thresh, y = 0, xend=top20_thresh, yend = max_value), colour = "black") +
    annotate("text", x = top20_thresh, y=max_value, label = paste0("top 20%\nOR: ", round(OR_top20$estimate[1], 1), sign_top20), vjust = -0.25) +
    geom_segment(aes(x=bottom20_thresh, y = 0, xend=bottom20_thresh, yend = max_value), colour = "black") +
    annotate("text", x = bottom20_thresh, y=max_value, label = paste0("bottom 20%\nOR: ", round(OR_bottom20$estimate[1], 1), sign_bot20), vjust = -0.25)
  Dens_plot
  
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_DENS_noComorbidity.svg"), Dens_plot, units = "cm", width = 13, height = 12, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_DENS_noComorbidity.pdf"), Dens_plot, units = "cm", width = 13, height = 12, dpi = 400)       
  
  # CALCULATE THE MEANS AND STANDARD ERRORS FOR EACH TIME POINT
  df_summary <- as.data.frame(DMR_categ %>%
                                group_by(Diagnosis) %>%
                                summarize(Mean = mean(AvgDMRmethylation, na.rm = TRUE),
                                          se = sd(AvgDMRmethylation, na.rm = TRUE) / sqrt(n()),
                                          median = median(AvgDMRmethylation, na.rm = TRUE),
                                          IQR = IQR(AvgDMRmethylation, na.rm = TRUE)))
  rownames(df_summary) <- levels(DMR_categ$Diagnosis)
  
  # CATEGORY PLOT
  Categ_plot <- ggplot(DMR_categ, aes(x = Diagnosis, y = AvgDMRmethylation, fill = Diagnosis)) +
    geom_violin(scale = "width", alpha = 0.3) + 
    geom_boxplot(width = 0.1, fill = "white") +
    geom_point(data = df_summary, aes(y = Mean, x = Diagnosis), color = "#14747b", size = 1) +
    th + th_transparent +
    ggtitle(GeneString) +
    labs(x = NULL, y = "Average DNAm of CpGs in DMR") + 
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a"))+
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
  Categ_plot
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_CATEG.pdf"), Categ_plot, units = "cm", width = 7, height = 10.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_CATEG.svg"), Categ_plot, units = "cm", width = 7, height = 10.5, dpi = 400)
  
  # CPG PLOTS
  
  cpg_df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 4))
  colnames(cpg_df) <- c("Basename", "Diagnosis", "cpg", "DNAm")
  
  DMR_cpgs_incl <- intersect(DMR_cpgs, rownames(beta_wobatch_noComorb))
  
  for (cg in DMR_cpgs_incl){
    single_cg_df <- data.frame(Basename = samples_noComorb$Basename, CaseCtrl = samples_noComorb$Diagnosis, cpg = cg, DNAm = NA)
    
    for (i in c(1:nrow(single_cg_df))){
      Basename <- single_cg_df$Basename[i]
      betas_cpg <- beta_wobatch_noComorb[cg, Basename]
      single_cg_df[i, "DNAm"] <- betas_cpg
    }
    
    cpg_df <- rbind(cpg_df, single_cg_df)
  }
  
  
  
  cpg_df$CaseCtrl <- factor(cpg_df$CaseCtrl, levels = c("CTRL", "OCD"))
  
  cpg_plot <- ggplot(cpg_df, aes(x=cpg, y = DNAm, fill = CaseCtrl)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    th + th_transparent + 
    ggtitle(GeneString) +
    labs(y = "DNA methylation") +
    scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(), 
          legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  cpg_plot
  
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_CPG.svg"), cpg_plot, units = "cm", width = 17, height = 13.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_CPG.pdf"), cpg_plot, units = "cm", width = 17, height = 13.5, dpi = 400)       
  
  
  
  median_data <- as.data.frame(cpg_df %>%
                                 group_by(cpg, CaseCtrl) %>%
                                 summarize(median = median(DNAm), .groups = "drop"))
  
  merged_cpg <- merge(x=median_data, y=anno, by.x = "cpg", by.y = "Name", all.x = TRUE)
  
  x_label <- paste0("position on ", chr)
  
  cpg_line_plot <- ggplot(merged_cpg, aes(x=MAPINFO, y = median, color = CaseCtrl, group = CaseCtrl)) +
    geom_point(size = 1) +
    geom_line() +
    labs(
      x = x_label,
      y = "median DNA methylation",
      color = NULL
    ) + th + th_transparent +
    ggtitle(GeneString) +
    scale_color_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
    theme(legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  cpg_line_plot
  
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_TRACK.svg"), cpg_line_plot, units = "cm", width = 20, height = 13.5, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_TRACK.pdf"), cpg_line_plot, units = "cm", width = 20, height = 13.5, dpi = 400)       
  
  
  median_diff <- unique(merged_cpg[, c("cpg", "CHR", "MAPINFO")])
  median_diff$median_OCD_minus_CTRL <- NA
  
  for (i in c(1:nrow(median_diff))){
    cpg <- median_diff$cpg[i]
    OCD_value <- merged_cpg$median[which(merged_cpg$cpg == cpg & merged_cpg$CaseCtrl == "OCD")]
    CTRL_value <- merged_cpg$median[which(merged_cpg$cpg == cpg & merged_cpg$CaseCtrl == "CTRL")]
    median_diff$median_OCD_minus_CTRL[i] <- OCD_value - CTRL_value
  }
  
  cpg_line_plot_diff <- ggplot(median_diff, aes(x=MAPINFO, y = median_OCD_minus_CTRL)) +
    geom_point(size = 1, color = "#891a1a") +
    geom_line(color = "#891a1a") +
    labs(
      x = x_label,
      y = "difference in median DNAm (OCD-CTRL)",
      color = NULL
    ) + th + th_transparent +
    ggtitle(GeneString) +
    theme(legend.position = "top",
          legend.box.background = element_blank(), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    geom_hline(color = "#084081", yintercept = 0)
  cpg_line_plot_diff
  
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_TRACK_DIFF.svg"), cpg_line_plot_diff, units = "cm", width = 20, height = 10, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/noComorbidity/DMR_CaseCtrl_", GeneString2, "_TRACK_DIFF.pdf"), cpg_line_plot_diff, units = "cm", width = 20, height = 10, dpi = 400)       
  
  
  
  
  
  
  
  
  
  Categ_plot <- Categ_plot + ggtitle("A")+ theme(plot.title = element_text(hjust = 0))
  Dens_plot <- Dens_plot + ggtitle("B")+ theme(plot.title = element_text(hjust = 0))
  #class(Dens_plot) <- class(Dens_plot)[class(Dens_plot) != "gg"]
  cpg_line_plot_diff <- cpg_line_plot_diff + ggtitle("C")+ theme(plot.title = element_text(hjust = 0))
  
  empty <- ggplot() + theme_void() #+ th + th_transparent
  
  design <- "
  ABBB
  CCCC
  DDDD
  "
  
  layout <- (Categ_plot + Dens_plot + empty + cpg_line_plot_diff) +
    plot_layout(design = design, heights = c(1, 0.1, 0.7)) +
    plot_annotation(title = GeneString,
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
  
  
  ggsave(paste0("output/Figures/DMRs/noComorbidity/combined/DMR_CaseCtrl_", GeneString2, "_COMB.svg"), layout, units = "cm", width = 22, height = 24, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/noComorbidity/combined/DMR_CaseCtrl_", GeneString2, "_COMB.pdf"), layout, units = "cm", width = 22, height = 24, dpi = 400)       
  
  rm(empty, layout)
  
}









################################################################################
# MQTL ANALYSIS
################################################################################

##### GENERAL THINGS ###########################################################

load("S:/Project/WP-genetics/05_Kira_mQTL/4_Final_SNP_files/Combined_samplesheet.RData")
load("S:/Project/WP-genetics/05_Kira_mQTL/4_Final_SNP_files/Combined_SNPmatrix.RData")
load("S:/Project/WP-genetics/05_Kira_mQTL/4_Final_SNP_files/Combined_map.RData")



stop("PERFORM PCA IN PLINK. This is not an error.")

snpsprobe_PCAscores <- read.table("S:/Project/WP-genetics/05_Kira_mQTL/5_PCA_results/CaseCtrl/PCA_result.eigenvec", header = FALSE)
colnames(snpsprobe_PCAscores) <- c("FID", "IID", paste0("PC", 1:20))

# filter based on successful genotyping PCs & order
samplesheet_comb_filt <- samplesheet_comb[which(samplesheet_comb$Resp_nr %in% snpsprobe_PCAscores$FID), ]
rownames(samplesheet_comb_filt) <- samplesheet_comb_filt$Resp_nr
samplesheet_comb_filt  <- samplesheet_comb_filt[sort(rownames(samplesheet_comb_filt)), ]

# extract genotyping PCs
rownames(snpsprobe_PCAscores) <- snpsprobe_PCAscores$FID
snpsprobe_PCAscores <- snpsprobe_PCAscores[sort(rownames(snpsprobe_PCAscores)), ]

snps_PC1 = snpsprobe_PCAscores$PC1
snps_PC2 = snpsprobe_PCAscores$PC2
snps_PC3 = snpsprobe_PCAscores$PC3
snps_PC4 = snpsprobe_PCAscores$PC4
snps_PC5 = snpsprobe_PCAscores$PC5
snps_PC6 = snpsprobe_PCAscores$PC6
snps_PC7 = snpsprobe_PCAscores$PC7
snps_PC8 = snpsprobe_PCAscores$PC8
snps_PC9 = snpsprobe_PCAscores$PC9
snps_PC10 = snpsprobe_PCAscores$PC10


# filter

for (i in c(1:nrow(samples))){
  Resp_nr <- samples$Resp_nr[i]
  if (Resp_nr == "no_Resp_nr"){
    new_resp_nr <- samples$Indiv_ID[i]
    samples$Resp_nr[i] <- new_resp_nr
  }
}

samples_filt <- samples[which(samples$Resp_nr %in% samplesheet_comb_filt$Resp_nr), ]
rownames(samples_filt) <- samples_filt$Resp_nr
samples_filt <- samples_filt[sort(rownames(samples_filt)), ]

beta_filt <- beta[, colnames(beta) %in% samples_filt$Basename]
beta_filt <- beta_filt[, samples_filt$Basename]
Mvalues_filt <- Mvalues[, colnames(Mvalues) %in% samples_filt$Basename]
Mvalues_filt <- Mvalues_filt[, samples_filt$Basename]

# covariables
age <- scale(as.numeric(samples_filt$Age))
smoking = scale(as.numeric(samples_filt$smoking))
sex <- as.factor(samples_filt$Sex)

# PCA
ctrl_sorted <- ctrl[, samples_filt$Basename]
ctrlprobe_PCAscores <- ctrl_probe_PCA(ctrl, samples_filt)
Ctrl_PC1 = ctrlprobe_PCAscores$PC1
Ctrl_PC2 = ctrlprobe_PCAscores$PC2
Ctrl_PC3 = ctrlprobe_PCAscores$PC3
Ctrl_PC4 = ctrlprobe_PCAscores$PC4
Ctrl_PC5 = ctrlprobe_PCAscores$PC5
Ctrl_PC6 = ctrlprobe_PCAscores$PC6
Ctrl_PC7 = ctrlprobe_PCAscores$PC7
Ctrl_PC8 = ctrlprobe_PCAscores$PC8
Ctrl_PC9 = ctrlprobe_PCAscores$PC9
Ctrl_PC10 = ctrlprobe_PCAscores$PC10
Ctrl_PC11 = ctrlprobe_PCAscores$PC11
Ctrl_PC12 = ctrlprobe_PCAscores$PC12
Ctrl_PC13 = ctrlprobe_PCAscores$PC13
Ctrl_PC14 = ctrlprobe_PCAscores$PC14
Ctrl_PC15 = ctrlprobe_PCAscores$PC15

cellcounts_PCAscores <- cell_type_PCA(samples_filt, c("Epi", "Fib", "comb_ICs"))
cell_PC1 = cellcounts_PCAscores$PC1
cell_PC2 = cellcounts_PCAscores$PC2


# ADAPT DNAm
# adjust beta values
covar_matrix <- as.matrix(data.frame(Age = age, 
                                     cell_PC1 = cell_PC1,
                                     cell_PC2 = cell_PC2, 
                                     smoking = smoking,
                                     Ctrl_PC1=Ctrl_PC1, 
                                     Ctrl_PC2=Ctrl_PC2, 
                                     Ctrl_PC3=Ctrl_PC3, 
                                     Ctrl_PC4=Ctrl_PC4, 
                                     Ctrl_PC5=Ctrl_PC5,
                                     Ctrl_PC6=Ctrl_PC6, 
                                     Ctrl_PC7=Ctrl_PC7, 
                                     Ctrl_PC8=Ctrl_PC8, 
                                     Ctrl_PC9=Ctrl_PC9, 
                                     Ctrl_PC10=Ctrl_PC10,
                                     Ctrl_PC11=Ctrl_PC11,
                                     Ctrl_PC12=Ctrl_PC12,
                                     Ctrl_PC13=Ctrl_PC13,
                                     Ctrl_PC14=Ctrl_PC14,
                                     Ctrl_PC15=Ctrl_PC15,
                                     snps_PC1=snps_PC1, 
                                     snps_PC2=snps_PC2,
                                     snps_PC3=snps_PC3,
                                     snps_PC4=snps_PC4,
                                     snps_PC5=snps_PC5,
                                     snps_PC6=snps_PC6,
                                     snps_PC7=snps_PC7,
                                     snps_PC8=snps_PC8,
                                     snps_PC9=snps_PC9,
                                     snps_PC10=snps_PC10))

Mvalues_wobatch_genotypes <- removeBatchEffect(Mvalues_filt, covariates = covar_matrix)
rm(covar_matrix)

beta_wobatch_genotypes <- M2B(Mvalues_wobatch_genotypes)
rm(Mvalues_wobatch_genotypes)


# MQTL ANALYSIS #

# NEW TABLE FOR RESULTS
mQTL_results_comb <- as.data.frame(matrix(data = NA, nrow =0, ncol = 10))
colnames(mQTL_results_comb) <- c("DMR", "SNP", 
                                 "Estimate_DNAm", "SE_DNAm", "t_value_DNAm", "p_value_DNAm",
                                 "Estimate_CaseCtrl", "SE_CaseCtrl", "z_value_CaseCtrl", "p_value_CaseCtrl")


for (j in c(1:nrow(dmrff_sign))){
  
  # EXTRACT INFOS FROM DMR
  chr <- dmrff_sign$chr[j]
  
  if (chr != "chrX" & chr != "chrY"){
    
    start <- as.numeric(dmrff_sign$start[j])
    end <- as.numeric(dmrff_sign$end[j])
    GeneString <- dmrff_sign$GeneString[j]
    GeneString2 <- gsub(";", "_", GeneString)
    DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
    beta_DMR <- beta_filt[rownames(beta_filt) %in% DMR_cpgs, ]
    
    # COMBINE DMR AVERAGE DNAm WITH SAMPLESHEET INFO
    DMR_contin <- add_avgDNAm(samples_filt, beta_DMR, "Basename", c("Resp_nr", "Diagnosis"))
    DMR_contin$Avg_Mvalues <-  B2M(DMR_contin[, "AvgDMRmethylation"])
    
    # SNPs FOR DMR
    SNPs <- Matt_map_overlapping$SNP[which(Matt_map_overlapping$chr == chr & 
                                             Matt_map_overlapping$pos >= (start - 50000) &
                                             Matt_map_overlapping$pos <= (end + 50000))]
    
    
    if(length(SNPs) > 0){
      
      # new results df
      mQTL_results <- as.data.frame(matrix(data = NA, nrow = length(SNPs), ncol = 10))
      rownames(mQTL_results) <- SNPs
      colnames(mQTL_results) <- c("DMR", "SNP", 
                                  "Estimate_DNAm", "SE_DNAm", "t_value_DNAm", "p_value_DNAm",
                                  "Estimate_CaseCtrl", "SE_CaseCtrl", "z_value_CaseCtrl", "p_value_CaseCtrl")
      mQTL_results$DMR <- GeneString
      
      
      
      # for each SNP
      for (k in c(1: length(SNPs))){
        SNP <- SNPs[k]
        SNP_genotypes <- as.vector(SNP_matrix_comb[DMR_contin$Resp_nr, SNP])
        DMR_contin$SNP <- SNP_genotypes
        DMR_contin$SNP <- as.numeric(DMR_contin$SNP)
        DMR_contin$SNP[which(DMR_contin$SNP == 0)] <- NA
        SNP_data <- DMR_contin$SNP
        SNP_data_frequency <- as.data.frame(table(SNP_data))
        CaseCtrl_status <- as.factor(DMR_contin$Diagnosis)
        
        
        # run linear model (no allele can be >95% of frequency!)
        if(nrow(SNP_data_frequency) > 1 & !any(SNP_data_frequency$Freq > 0.95*sum(SNP_data_frequency$Freq))){
          
          mQTL_results$SNP[k] <- SNP
          
          # DNAm ~ genotype
          lm_result <- summary(lm(DMR_contin$Avg_Mvalues ~ SNP_data + sex + age + smoking + 
                                    cell_PC1 + cell_PC2 + 
                                    Ctrl_PC1 + Ctrl_PC2 + Ctrl_PC3 + Ctrl_PC5 + Ctrl_PC6 + Ctrl_PC7 + Ctrl_PC8 + Ctrl_PC9 + Ctrl_PC10 +
                                    Ctrl_PC11 + Ctrl_PC12 + Ctrl_PC13 + Ctrl_PC14 + Ctrl_PC15 +
                                    snps_PC1 + snps_PC2 + snps_PC3 + snps_PC4 + snps_PC5 + snps_PC6 + snps_PC7 + snps_PC8 + snps_PC9  + snps_PC10))
          
          glm_result_CC <- summary(glm(CaseCtrl_status ~ SNP_data + sex + age +
                                    snps_PC1 + snps_PC2 + snps_PC3, family = "binomial"))
          
          
          if ("SNP_data" %in% rownames(lm_result$coefficients)){
            mQTL_results$Estimate_DNAm[k] <- lm_result$coefficients["SNP_data", "Estimate"]
            mQTL_results$SE_DNAm[k] <- lm_result$coefficients["SNP_data", "Std. Error"]
            mQTL_results$t_value_DNAm[k] <- lm_result$coefficients["SNP_data", "t value"]
            mQTL_results$p_value_DNAm[k] <- lm_result$coefficients["SNP_data", "Pr(>|t|)"]
            
            mQTL_results$Estimate_CaseCtrl[k] <- glm_result_CC$coefficients["SNP_data", "Estimate"]
            mQTL_results$SE_CaseCtrl[k] <- glm_result_CC$coefficients["SNP_data", "Std. Error"]
            mQTL_results$z_value_CaseCtrl[k] <- glm_result_CC$coefficients["SNP_data", "z value"]
            mQTL_results$p_value_CaseCtrl[k] <- glm_result_CC$coefficients["SNP_data", "Pr(>|z|)"]
          }
        }
        rm(SNP_data)
      }
      mQTL_results_comb <- rbind(mQTL_results_comb, mQTL_results)
    }
  }
}


# CORRECTION FOR MULTIPLE TESTING
mQTL_results_comb <- mQTL_results_comb[which(!is.na(mQTL_results_comb$SNP)), ]
mQTL_results_comb$Bonferroni_all_DNAm <- p.adjust(p = mQTL_results_comb$p_value_DNAm, method = "bonferroni")

# EXTRACT SIGNIFICANT mQTLS AND ADJUST P VALUES
mQTL_results_comb_sign <- mQTL_results_comb[which(mQTL_results_comb$Bonferroni_all_DNAm < 0.05), ]

rm(mQTL_results)


# ADJUST P VALUE
write_xlsx(mQTL_results_comb_sign, path = "output/Tables/mQTL_sign.xlsx")
save(mQTL_results_comb_sign, file = "output/RData/mQTL_sign.RData")

write_xlsx(mQTL_results_comb, path = "output/Tables/mQTL.xlsx")
save(mQTL_results_comb, file = "output/RData/mQTL.RData")



# Clumping

summary_stats_plink <- data.frame(SNP = rownames(mQTL_results_comb_sign), CHR = NA, BP = NA, P = mQTL_results_comb_sign$p_value_CaseCtrl)

for (i in c(1:nrow(summary_stats_plink))){
  SNP <- summary_stats_plink$SNP[i]
  chr <- as.numeric(Matt_map_overlapping$chr_numeric[which(Matt_map_overlapping$SNP == SNP)])
  pos <- as.numeric(Matt_map_overlapping$pos[which(Matt_map_overlapping$SNP == SNP)])
  summary_stats_plink$CHR[i] <- chr
  summary_stats_plink$BP[i] <- pos
}


write.table(summary_stats_plink, "S:/Project/WP-genetics/05_Kira_mQTL/5_Clumping/CaseCtrl/association_results.assoc",
            quote = FALSE, row.names = FALSE, sep = "\t")

stop("Clumping - not an error")
# 30 clumps formed from 491 top variants

clumped_snps <- read.table("S:/Project/WP-genetics/05_Kira_mQTL/5_Clumping/CaseCtrl/clumped_snps_final.clumped", 
                           header = TRUE, stringsAsFactors = FALSE)

mQTL_results_comb_clumped <- mQTL_results_comb_sign[which(mQTL_results_comb_sign$SNP %in% clumped_snps$SNP), ]

mQTL_results_comb_clumped$adj_p_CaseCtrl_Bonferroni <- p.adjust(mQTL_results_comb_clumped$p_value_CaseCtrl, method = "bonferroni")

write_xlsx(mQTL_results_comb_clumped, path = "output/Tables/mQTL_sign_clumped.xlsx")
save(mQTL_results_comb_clumped, file = "output/RData/mQTL_sign_clumped.RData")


# VISUALISATION
for (j in c(1:nrow(mQTL_results_comb_clumped))){
  gene <- mQTL_results_comb_clumped$DMR[j]
  snp <- mQTL_results_comb_clumped$SNP[j]
  
  # EXTRACT INFOS FROM DMR
  chr <- dmrff_sign$chr[which(dmrff_sign$GeneString == gene)]
  start <- dmrff_sign$start[which(dmrff_sign$GeneString == gene)]
  end <- dmrff_sign$end[which(dmrff_sign$GeneString == gene)]
  GeneString <- dmrff_sign$GeneString[which(dmrff_sign$GeneString == gene)]
  GeneString2 <- gsub(";", "_", GeneString)
  
  # EXTRACT P VALUES
  snp_pvalue_DNAm <- mQTL_results_comb_clumped$BH_all_DNAm[which(mQTL_results_comb_clumped$SNP == snp)]
  snp_pvalue_CaseCtrl <- mQTL_results_comb_clumped$adj_p_CaseCtrl_Bonferroni[which(mQTL_results_comb_clumped$SNP == snp)]
  
  
  ####### DNA METHYLATION PLOT ########
  
  # EXTRACT ADAPTED BETA VALUES
  DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
  beta_DMR <- beta_wobatch_genotypes[rownames(beta_wobatch_genotypes) %in% DMR_cpgs, ]
  
  # COMBINE DMR AVERAGE DNAm WITH SAMPLESHEET INFO
  DMR_contin <- add_avgDNAm(samples_filt, beta_DMR, "Basename", c("Resp_nr", "Diagnosis"))
  
  DMR_contin$Avg_Mvalues <-  B2M(DMR_contin[, "AvgDMRmethylation"])
  
  SNP_genotypes <- SNP_matrix_comb_filt[DMR_contin$Resp_nr, snp]
  
  #combine SNP data with DNAm data
  comb_DMR_SNP <- cbind(DMR_contin, SNP_genotypes)
  comb_DMR_SNP[, 6] <- as.numeric(comb_DMR_SNP[, 6])
  comb_DMR_SNP[, 6][comb_DMR_SNP[, 6] == 0] <- NA
  comb_DMR_SNP <- comb_DMR_SNP[, c(1, 6)]
  colnames(comb_DMR_SNP) <- c("Basename", "SNP1")
  comb_DMR_SNP <- merge(comb_DMR_SNP, DMR_contin, by = "Basename", all.x = TRUE)
  
  nr_genotypes <- length(unique(comb_DMR_SNP$SNP1[!is.na(comb_DMR_SNP$SNP1)]))
  
  # add genotypes
  alleles_SNP1 <- Matt_map_overlapping[Matt_map_overlapping$SNP == snp, c("allele1", "allele2")]
  sorted_alleles_SNP1 <- sort(unlist(as.vector(alleles_SNP1[1, ])))
  comb_DMR_SNP$SNP1[which(comb_DMR_SNP$SNP1 == 1)] <- paste0(alleles_SNP1$allele1[1], alleles_SNP1$allele1[1])
  comb_DMR_SNP$SNP1[which(comb_DMR_SNP$SNP1 == 2)] <- paste0(sorted_alleles_SNP1[1], sorted_alleles_SNP1[2])
  comb_DMR_SNP$SNP1[which(comb_DMR_SNP$SNP1 == 3)] <- paste0(alleles_SNP1$allele2[1], alleles_SNP1$allele2[1])
  
  # remove missing entries
  comb_DMR_SNP_woNA <- comb_DMR_SNP[which(!is.na(comb_DMR_SNP$SNP1)), ]
  comb_DMR_SNP_woNA <- comb_DMR_SNP_woNA[order(comb_DMR_SNP_woNA$SNP1), ]
  rm(comb_DMR_SNP)
  
  # order the levels of genotypes
  right_order <- c("AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT")
  level_order <- intersect(right_order, comb_DMR_SNP_woNA$SNP1)
  comb_DMR_SNP_woNA$SNP1 <- factor(comb_DMR_SNP_woNA$SNP1, levels = level_order)
  
  rm(alleles_SNP1, sorted_alleles_SNP1)
  
  
  #stars annotation p values
  if (snp_pvalue_DNAm < 0.05){
    text_DNAm <- "*"
  } 
  if (snp_pvalue_DNAm < 0.01) {
    text_DNAm <- "**"
  }
  if (snp_pvalue_DNAm < 0.001) {
    text_DNAm <- "***"
  } 
  if (snp_pvalue_DNAm < 0.0001) {
    text_DNAm <- "****"
  } 
  if (snp_pvalue_DNAm >= 0.05) {
    text_DNAm <- "ns"
  } 
  
  
  
  
  diff <- max(comb_DMR_SNP_woNA$AvgDMRmethylation) - min(comb_DMR_SNP_woNA$AvgDMRmethylation)
  
  if (nr_genotypes == 3){
    # DNA methylation plot
    DNAm_SNP_plot <- ggplot(comb_DMR_SNP_woNA, aes(x=SNP1, y = AvgDMRmethylation, fill = SNP1)) +
      geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8), fill = "lightgrey") +
      th + th_transparent +
      ggtitle(paste0(GeneString, "\n", snp)) +
      labs(y = "average CpG DNAm in DMR", x = "genotype") +
      theme(legend.position = "none",
            legend.box.background = element_blank(), 
            strip.text = element_text(face = "bold", size = 10), 
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.spacing = unit(1.5, "lines")) +
      annotate("segment", x = 1, xend = 3, y = max(comb_DMR_SNP_woNA$AvgDMRmethylation) + diff * 0.09, 
               yend = max(comb_DMR_SNP_woNA$AvgDMRmethylation) + diff * 0.09, color = "black") +
      annotate("text", x = 2, y = max(comb_DMR_SNP_woNA$AvgDMRmethylation) + diff * 0.13, label = text_DNAm, size = 5)
  }
  
  if (nr_genotypes == 2){
    DNAm_SNP_plot <- ggplot(comb_DMR_SNP_woNA, aes(x=SNP1, y = AvgDMRmethylation, fill = SNP1)) +
      geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8), fill = "lightgrey") +
      th + th_transparent +
      ggtitle(paste0(GeneString, "\n", snp)) +
      labs(y = "average CpG DNAm in DMR", x = "genotype") +
      theme(legend.position = "none",
            legend.box.background = element_blank(), 
            strip.text = element_text(face = "bold", size = 10), 
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.spacing = unit(1.5, "lines")) +
      annotate("segment", x = 1, xend = 2, y = max(comb_DMR_SNP_woNA$AvgDMRmethylation) + diff * 0.09, 
               yend = max(comb_DMR_SNP_woNA$AvgDMRmethylation) + diff * 0.09, color = "black") +
      annotate("text", x = 1.5, y = max(comb_DMR_SNP_woNA$AvgDMRmethylation) + diff * 0.115, label = text_DNAm, size = 5)
  }  
  
  
  # P value stars
  if (snp_pvalue_CaseCtrl < 0.05){
    text_CaseCtrl <- "*"
  } 
  if (snp_pvalue_CaseCtrl < 0.01) {
    text_CaseCtrl <- "**"
  }
  if (snp_pvalue_CaseCtrl < 0.001) {
    text_CaseCtrl <- "***"
  } 
  if (snp_pvalue_CaseCtrl < 0.0001) {
    text_CaseCtrl <- "****"
  }
  if (snp_pvalue_CaseCtrl >= 0.05) {
    text_CaseCtrl <- "ns"
  }
  
  
  df_percentage <- as.data.frame(comb_DMR_SNP_woNA %>%
    group_by(SNP1, Diagnosis) %>%
    summarise(Count = n(), .groups = "drop"))
  
  df_percentage <- as.data.frame(df_percentage %>%
                                   group_by(SNP1) %>%
                                   mutate(Percentage = Count/ sum(Count) * 100))
  
  # sum_CTRL <- sum(df_percentage$Count[which(df_percentage$Diagnosis == "CTRL")])
  # sum_OCD <- sum(df_percentage$Count[which(df_percentage$Diagnosis == "OCD")])
  # 
  # df_percentage$percentage <- NA
  # 
  # for (i in 1:nrow(df_percentage)){
  #   casectrl <- df_percentage$Diagnosis[i]
  #   if (casectrl == "CTRL"){
  #     df_percentage$percentage[i] <- (df_percentage$Count[i] / sum_CTRL) *100
  #   } 
  #   if (casectrl == "OCD"){
  #     df_percentage$percentage[i] <- (df_percentage$Count[i] / sum_OCD) *100
  #   } 
  # }
  
  if (nr_genotypes == 3){
    
    Diagnosis_SNP_plot <- ggplot(df_percentage, aes(x = SNP1, y = Percentage, fill = Diagnosis)) + 
      geom_bar(stat = "identity", position = "stack", alpha = 0.8, color = "black", linewidth = 0.5) + 
      labs(title = paste0(GeneString2, "\n",snp), x = "genotype", y = "% OCD/CTRL per genotype") +
      th + th_transparent +
      scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
      theme(legend.key = element_blank(),
            legend.background = element_blank(),
            legend.position = "top",
            legend.box.background = element_blank(), 
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, 112), breaks = c(0,20,40,60,80,100)) +
      annotate("segment", x = 1, xend = 3, y = 103, yend = 103, color = "black") +
      annotate("text", x = 2, y = 108, label = text_CaseCtrl, size = 5)
    Diagnosis_SNP_plot
  }
  
  if (nr_genotypes == 2){
    
    Diagnosis_SNP_plot <- ggplot(df_percentage, aes(x = SNP1, y = Percentage, fill = Diagnosis)) + 
      geom_bar(stat = "identity",position = "stack", alpha = 0.8, color = "black", linewidth = 0.5) + 
      labs(title = paste0(GeneString2, "\n",snp), x = "genotype", y = "% OCD/CTRL per genotype") +
      th + th_transparent +
      scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
      theme(legend.key = element_blank(),
            legend.background = element_blank(),
            legend.position = "top",
            legend.box.background = element_blank(), 
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, 112), breaks = c(0,20,40,60,80,100)) +
      annotate("segment", x = 1, xend = 2, y = 103, yend = 103, color = "black") +
      annotate("text", x = 1.5, y =108, label = text_CaseCtrl, size = 5)
    Diagnosis_SNP_plot
  }
  
  
  if (nr_genotypes == 3){
    
    Diagnosis_SNP_count_plot <- ggplot(df_percentage, aes(x = SNP1, y = Count, fill = Diagnosis)) + 
      geom_bar(stat = "identity", position = position_dodge(preserve = "single"), alpha = 0.8, color = "black", linewidth = 0.5) + 
      labs(title = paste0(GeneString2, "\n",snp), x = "genotype", y = "OCD/CTRL count") +
      th + th_transparent +
      scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
      theme(legend.key = element_blank(),
            legend.background = element_blank(),
            legend.position = "top",
            legend.box.background = element_blank(), 
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, (max(df_percentage$Count) + 20))) +
      annotate("segment", x = 1, xend = 3, y = max(df_percentage$Count) + 5, 
               yend = max(df_percentage$Count) + 5, color = "black") +
      annotate("text", x = 2, y = max(df_percentage$Count) + 15, label = text_CaseCtrl, size = 5)
    Diagnosis_SNP_count_plot
  }
  
  if (nr_genotypes == 2){
    
    Diagnosis_SNP_count_plot <- ggplot(df_percentage, aes(x = SNP1, y = Percentage, fill = Diagnosis)) + 
      geom_bar(stat = "identity", position = position_dodge(preserve = "single"), alpha = 0.8, color = "black", linewidth = 0.5) + 
      labs(title = paste0(GeneString2, "\n",snp), x = "genotype", y = "OCD/CTRL count") +
      th + th_transparent +
      scale_fill_manual(values = c("CTRL" = "#084081", "OCD" = "#891a1a")) +
      theme(legend.key = element_blank(),
            legend.background = element_blank(),
            legend.position = "top",
            legend.box.background = element_blank(), 
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, (max(df_percentage$Percentage) + 20))) +
      annotate("segment", x = 1, xend = 2, y = max(df_percentage$Percentage) + 5, 
               yend = max(df_percentage$Percentage) + 5, color = "black") +
      annotate("text", x = 1.5, y = max(df_percentage$Percentage) + 15, label = text_CaseCtrl, size = 5)
    Diagnosis_SNP_count_plot
  }
  
  
  #save
  ggsave(paste0("output/Figures/DMRs/mQTL/DMR_CaseCtrl_", GeneString2, "_", snp, "_DNAm_SNP.svg"), DNAm_SNP_plot, units = "cm", width = 6.5, height = 11, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/mQTL/DMR_CaseCtrl_", GeneString2, "_", snp, "_DNAm_SNP.pdf"), DNAm_SNP_plot, units = "cm", width = 6.5, height = 11, dpi = 400)       
  save(DNAm_SNP_plot, diff, file = paste0("output/Figures/DMRs/mQTL/DMR_CaseCtrl_", GeneString2, "_", snp, "_DNAm_SNP.RData"))
  ggsave(paste0("output/Figures/DMRs/mQTL/DMR_CaseCtrl_", GeneString2, "_", snp, ".svg"), Diagnosis_SNP_plot, units = "cm", width = 6.5, height = 11, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/mQTL/DMR_CaseCtrl_", GeneString2, "_", snp, ".pdf"), Diagnosis_SNP_plot, units = "cm", width = 6.5, height = 11, dpi = 400)       
  save(Diagnosis_SNP_plot, df_percentage, file = paste0("output/Figures/DMRs/mQTL/DMR_CaseCtrl_", GeneString2, "_", snp, ".RData"))
  ggsave(paste0("output/Figures/DMRs/mQTL/DMR_CaseCtrl_", GeneString2, "_", snp, "_count.svg"), Diagnosis_SNP_count_plot, units = "cm", width = 6.5, height = 11, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/mQTL/DMR_CaseCtrl_", GeneString2, "_", snp, "_count.pdf"), Diagnosis_SNP_count_plot, units = "cm", width = 6.5, height = 11, dpi = 400)       
  save(Diagnosis_SNP_count_plot, df_percentage, file = paste0("output/Figures/DMRs/mQTL/DMR_CaseCtrl_", GeneString2, "_", snp, "_count.RData"))
  
  
  # combined plot
  DNAm_SNP_plot <- DNAm_SNP_plot + ggtitle("A") + theme(plot.title = element_text(hjust = 0))
  Diagnosis_SNP_count_plot <- Diagnosis_SNP_count_plot + ggtitle("B")+ theme(plot.title = element_text(hjust = 0))
  Diagnosis_SNP_plot <- Diagnosis_SNP_plot + ggtitle("C")+ theme(plot.title = element_text(hjust = 0))
  empty <- ggplot() + theme_void()
  empty2 <- ggplot() + theme_void()
  
  design <- "AAAAABCCCCCCCDEEEEEE"
  
  layout <- (DNAm_SNP_plot + empty + Diagnosis_SNP_count_plot + empty2 + Diagnosis_SNP_plot) +
    plot_layout(design = design) +
    plot_annotation(title = paste0(GeneString, "\n", snp),
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
  
  
  ggsave(paste0("output/Figures/DMRs/mQTL/combined/DMR_CaseCtrl_", GeneString2, "_", snp, "_COMB.svg"), layout, units = "cm", width = 19, height = 12, dpi = 400)
  ggsave(paste0("output/Figures/DMRs/mQTL/combined/DMR_CaseCtrl_", GeneString2, "_", snp, "_COMB.pdf"), layout, units = "cm", width = 19, height = 12, dpi = 400)       
  
}


table(mQTL_results_comb_clumped$DMR)

# ARHGAP45           ARHGEF17              BRSK2              CCDC8               CDX1             CYB5R3              DLEU7               FAT3 
# 2                  1                  4                  8                  2                  7                  3                  3 
# FXYD1 JAKMIP1-DT;JAKMIP1              KIFC3          LINC02580               NOM1               PAX9              PGBD5              RGS12 
# 1                  1                  1                  5                  4                  1                  7                  6 
# RUNX3       UNANNOTATED3            ZNF833P 
# 1                  3                  4 

