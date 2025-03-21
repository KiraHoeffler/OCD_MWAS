
################################################################################
# SETUP
################################################################################

# SET PATHS
dir_gen <- "S://Project/WP-epigenetics/10_CaseCtrl_EWAS/" #working directory
ctrl_path <- "S://Project/WP-epigenetics/04_Pipeline_PhaseI/output/RData/Positive_ctrlprobe_intensities.RData"
SNP_path <- "S://Project/WP-epigenetics/04_Pipeline_PhaseI/output/RData/comb_SNPs.RData"

# SET COLOURS
epi_color <- "#489d5e"
fib_color <- "#73bf86"
B_color <-  "#a8ddb5"
NK_color <- "#7bccc4"
CD4T_color <- "#4eb3d3"
CD8T_color <- "#2b8cbe"
mono_color <- "#0868ac"
eosino_color <- "#084081"
neutro_color <- "#052448"

ctrl_color <- "#084081"
case_color <- "#891a1a"


################################################################################
# LOAD LIBRARIES
################################################################################

# LOAD LIBRARIES
library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
library(shiny)
library(readxl)
library(ggplot2)
library(writexl)
library(reshape2)
library(dplyr)
library(patchwork)


################################################################################
# IMPORT
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

#New Folders
if (!dir.exists("output/Tables/CelltypeProportions")) dir.create("output/Tables/CelltypeProportions")
if (!dir.exists("output/Figures/CelltypeProportions")) dir.create("output/Figures/CelltypeProportions")


# LOAD SAMPLESHEET 
samples <- as.data.frame(read_excel("output/Tables/Final_samplesheet_CaseCtrl.xlsx"))

# LOAD THEME
load("output/RData/theme.RData")
load("output/RData/theme_transparent.Rdata")

# CTRL PROBES
load(ctrl_path) 

# SNP PROBES
load(SNP_path) 

# FUNCTIONS
source("../05_Analysis_old/0_General/Functions/PCAs_for_model.R")
source("../05_Analysis_old/0_General/Functions/cell_prop_analyses.R")

################################################################################
# GENERAL
################################################################################

# CELL TYPES INCLUDED IN THE ANALYSES (THAT ARE NOT ZERO)
cell_type_vector <- c("Epi", "Fib", "B", "NK", "CD4T", "Mono", "Neutro")
cell_type_labels <- c("Epithelial cells", "Fibroblasts", "B cells", "NK cells", "CD4T cells", "Monocytes", "Neutrophils")

################################################################################
# Baseline
################################################################################

# PCAs FOR MODELS
ctrl_PCA <- ctrl_probe_PCA(ctrl, samples)
SNP_PCA <- ancestry_PCA(SNP_info = comb_SNPs, samplesheet = samples)

# EXTRACT WHAT IS NEEDED TO RUN MODEL
variables_df <- data.frame(
  Diagnosis = as.factor(samples$Diagnosis),
  Sex  = as.factor(samples$Sex),
  Age = scale(as.numeric(samples$Age)),
  smoking = scale(as.numeric(samples$smoking)),
  Epi = scale(as.numeric(samples$Epi)),
  Fib = scale(as.numeric(samples$Fib)),
  B = scale(as.numeric(samples$B)),
  NK = scale(as.numeric(samples$NK)),
  CD4T = scale(as.numeric(samples$CD4T)),
  Mono = scale(as.numeric(samples$Mono)),
  Neutro = scale(as.numeric(samples$Neutro)),
  ctrl_PC1 = ctrl_PCA$PC1,
  ctrl_PC2 = ctrl_PCA$PC2,
  ctrl_PC3 = ctrl_PCA$PC3,
  ctrl_PC5 = ctrl_PCA$PC5,
  ctrl_PC6 = ctrl_PCA$PC6,
  ctrl_PC7 = ctrl_PCA$PC7,
  ctrl_PC8 = ctrl_PCA$PC8,
  ctrl_PC9 = ctrl_PCA$PC9,
  ctrl_PC10 = ctrl_PCA$PC10,
  ctrl_PC11 = ctrl_PCA$PC11,
  ctrl_PC12 = ctrl_PCA$PC12,
  ctrl_PC13 = ctrl_PCA$PC13,
  ctrl_PC14 = ctrl_PCA$PC14,
  ctrl_PC15 = ctrl_PCA$PC15,
  snps_PC1 = SNP_PCA$PC1,
  snps_PC2 = SNP_PCA$PC2
  )

# RUN MODEL AND SAVE RESULTS

# combine cell types in list
cell_type_list <- list()
for (cellprop in cell_type_vector){
  cell_type_list[[cellprop]] <- variables_df[, cellprop]
}

# create empty df for results
comb_results <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 17))
colnames(comb_results) <- c("Cell_type", 
                            "Estimate_Diagnosis", "Std_Error_Diagnosis", "t_value_Diagnosis", "p_value_Diagnosis",
                            "Estimate_Sex", "Std_Error_Sex", "t_value_Sex", "p_value_Sex",
                            "Estimate_Age", "Std_Error_Age", "t_value_Age", "p_value_Age", 
                            "Estimate_smoking", "Std_Error_smoking", "t_value_smoking", "p_value_smoking")

# run model and add results
for (i in c(1:length(cell_type_vector))){
  cell_type <- cell_type_list[[i]]
  
  variables_df_copy <- variables_df
  variables_df_copy$cell_type <- cell_type
  
  model <- lm(cell_type ~ Diagnosis + Sex + Age + smoking + 
                ctrl_PC1 + ctrl_PC2 + ctrl_PC3 + ctrl_PC5 + ctrl_PC6 + ctrl_PC7 + ctrl_PC8 + ctrl_PC9 + ctrl_PC10 +
                ctrl_PC11 + ctrl_PC12 + ctrl_PC13 + ctrl_PC14 + ctrl_PC15 +
                snps_PC1 + snps_PC2, data = variables_df_copy)
  model_summary <- summary(model)
  rm(variables_df_copy)
  
  results <- data.frame(Cell_type = cell_type_vector[i],
                        Estimate_Diagnosis = model_summary$coefficients["DiagnosisOCD", "Estimate"],
                        Std_Error_Diagnosis = model_summary$coefficients["DiagnosisOCD", "Std. Error"],
                        t_value_Diagnosis = model_summary$coefficients["DiagnosisOCD", "t value"],
                        p_value_Diagnosis = model_summary$coefficients["DiagnosisOCD", "Pr(>|t|)"],
                        Estimate_Sex = model_summary$coefficients["SexM", "Estimate"],
                        Std_Error_Sex = model_summary$coefficients["SexM", "Std. Error"],
                        t_value_Sex = model_summary$coefficients["SexM", "t value"],
                        p_value_Sex = model_summary$coefficients["SexM", "Pr(>|t|)"],
                        Estimate_Age = model_summary$coefficients["Age", "Estimate"],
                        Std_Error_Age = model_summary$coefficients["Age", "Std. Error"],
                        t_value_Age = model_summary$coefficients["Age", "t value"],
                        p_value_Age = model_summary$coefficients["Age", "Pr(>|t|)"],
                        Estimate_smoking =  model_summary$coefficients["smoking", "Estimate"],
                        Std_Error_smoking = model_summary$coefficients["smoking", "Std. Error"],
                        t_value_smoking = model_summary$coefficients["smoking", "t value"],
                        p_value_smoking = model_summary$coefficients["smoking", "Pr(>|t|)"])
  
  comb_results <- rbind(comb_results, results)
  
}

################################################################################
# Number of independent tests (Li/Ji method)
################################################################################

samples_red <- samples[, cell_type_vector]
cor_matrix <- cor(samples_red, use = "pairwise.complete.obs")
eigen_values <- eigen(cor_matrix, symmetric = TRUE)$values
M_eff <- sum(eigen_values/ (eigen_values + 1))
print(M_eff) #1.82269

comb_results$adj_p_Diagnosis <- comb_results$p_value_Diagnosis * M_eff
comb_results$adj_p_Sex <- comb_results$p_value_Sex * M_eff
comb_results$adj_p_Age <- comb_results$p_value_Age * M_eff
comb_results$adj_p_Smoking <- comb_results$p_value_smoking * M_eff

write_xlsx(comb_results, path = "output/Tables/CelltypeProportions/CellTypes_CaseCtrl.xlsx")






################################################################################
# VISUALISATION
################################################################################

################### STACK PLOT #################################################

# MAKE OVERVIEW TABLE
cell_type_prop_hepi <- as.data.frame(samples[, c("Epi", "Fib", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino")])

se <- function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))) # define function standard error
summary_table <- data.frame(
  mean = sapply(cell_type_prop_hepi, mean),
  SD = sapply(cell_type_prop_hepi, sd),
  SE = sapply(cell_type_prop_hepi, se)
)
summary_table$Cell_type <- rownames(summary_table)

#EXPORT
write_xlsx(summary_table, path = "output/Tables/CelltypeProportions/Cell_type_proportions.xlsx")

# CREATE LONG DATAFRAME
cell_type_prop_hepi$SampleName <- samples$Basename
cell_type_prop_hepi_long <- melt(cell_type_prop_hepi, id.vars = "SampleName",
                                 variable.name = "CellType", value.name = "EstimatedProportion")
cell_type_prop_hepi_long$CellType <- as.character(cell_type_prop_hepi_long$CellType)
cell_type_prop_hepi_long$CellType <- factor(cell_type_prop_hepi_long$CellType, levels = c("Epi", "Fib", "B", "NK", "CD4T", "CD8T", "Mono", "Eosino", "Neutro"))

# ORDERING SAMPLES BASED ON EPITHELIAL CELL TYPE PROPORTIONS
ordering <- cell_type_prop_hepi_long %>%
  filter(CellType == "Epi") %>%
  arrange(desc(EstimatedProportion)) %>%
  .$SampleName
cell_type_prop_hepi_long$SampleName <- factor(cell_type_prop_hepi_long$SampleName, levels = ordering)

# PLOT CELL TYPE PROPORTIONS
celltype_prop_hepi_plot <- ggplot(cell_type_prop_hepi_long, aes(x = SampleName, y = EstimatedProportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) + th + th_transparent + 
  labs(y = "Estimated cell type proportion", fill = "Cell type") +
  xlab(NULL) + 
  theme(axis.text.x = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        axis.line.x = element_blank(), 
        legend.box.background = element_blank(), 
        legend.key = element_blank()) +
  scale_fill_manual(values = c("Epi" = epi_color, "Fib" = fib_color, "B" = B_color, "NK" = NK_color, 
                               "CD4T" = CD4T_color, "CD8T" = CD8T_color, "Mono" = mono_color, 
                               "Eosino" = eosino_color, "Neutro" = neutro_color))
celltype_prop_hepi_plot

# EXPORT PLOT
ggsave("output/Figures/CelltypeProportions/Cell_type_proportions_hepi.pdf", celltype_prop_hepi_plot, units = "cm", width = 11, height = 8, dpi = 400)
ggsave("output/Figures/CelltypeProportions/Cell_type_proportions_hepi.svg", celltype_prop_hepi_plot, units = "cm", device = "svg",  width = 11, height = 8, dpi = 400)
save(celltype_prop_hepi_plot, file = "output/Figures/CelltypeProportions/Cell_type_proportions_hepi.RData")

########################## VIOLIN PLOT #########################################

# SORT LONG DATA LEVELS OF CELL TYPE AFTER MEAN OF CELL TYPE PROPORTIONS
cell_type_prop_hepi2 <- cell_type_prop_hepi[, -ncol(cell_type_prop_hepi)]
means_hepi <- data.frame(colMeans(cell_type_prop_hepi2))
colnames(means_hepi) <- "mean_values"
means_hepi$CellTypes <- rownames(means_hepi)
means_hepi_sort <- means_hepi[order(means_hepi$mean_values, decreasing = TRUE), ]
cell_type_prop_hepi_long2 <- cell_type_prop_hepi_long
cell_type_prop_hepi_long2$CellType <- factor(cell_type_prop_hepi_long2$CellType, levels = means_hepi_sort$CellTypes)

# VIOLIN PLOT
Hepi_violin_plot <- ggplot(cell_type_prop_hepi_long2, aes(x = CellType, y = EstimatedProportion, fill = CellType)) +
  geom_violin(scale = "width", alpha = 0.8) + th + th_transparent + 
  geom_boxplot(width = 0.1, fill = "white") + 
  labs(x = "Cell type", y = "Estimated cell type proportion", fill = "Cell type") +
  scale_fill_manual(values = c("Epi" = epi_color, "Fib" = fib_color, "B" = B_color, "NK" = NK_color, 
                               "CD4T" = CD4T_color, "CD8T" = CD8T_color, "Mono" = mono_color, 
                               "Eosino" = eosino_color, "Neutro" = neutro_color)) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1.05))
Hepi_violin_plot

ggsave("output/Figures/CelltypeProportions/Cell_type_proportions_hepi_violin.pdf", Hepi_violin_plot, units = "cm", width = 15, height = 12, dpi = 400)
ggsave("output/Figures/CelltypeProportions/Cell_type_proportions_hepi_violin.svg", Hepi_violin_plot, units = "cm", device = "svg", width = 15, height = 12, dpi = 400)
save(Hepi_violin_plot, file = "output/Figures/CelltypeProportions/Cell_type_proportions_hepi_violin.RData")



################## CATEG_CELL_PLOT ##############################################



categ_cell_plot <- function(samples, cell_type, cell_type_name){
  samples_red <- samples[, c("Diagnosis", cell_type)]
  colnames(samples_red) <- c("Diagnosis", "specific_cellprop")
  
  Violin_plot <- ggplot(samples_red, aes(x = Diagnosis, y = specific_cellprop, fill = Diagnosis)) +
    geom_violin(scale = "width", alpha = 0.5) + th + th_transparent + 
    geom_boxplot(width = 0.1, fill = "white") + 
    labs(y = "Estimated cell type proportion", x = NULL) +
    ggtitle(cell_type_name) +
    scale_fill_manual(values = c("CTRL" = ctrl_color, "OCD" = case_color)) +
    theme(axis.title = element_text(face = "bold"),
          axis.title.y = element_text(margin = margin(r = 10)),
          legend.text = element_text(size = 12),
          legend.title = element_text(face = "bold"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  return(Violin_plot)
}


for (i in c(1:length(cell_type_vector))){
  cell_type <- cell_type_vector[i]
  cell_type_name <- cell_type_labels[i]
  
  plot <- categ_cell_plot(samples = samples, cell_type = cell_type, cell_type_name = cell_type_name)
  
  pvalue <- comb_results$adj_p_Diagnosis[which(comb_results$Cell_type == cell_type)]
  
  if(pvalue > 0.05){
    star <- "ns"
  }
  if(pvalue <= 0.05){
    star <- "*"
  }
  if(pvalue <= 0.01){
    star <- "**"
  }
  if(pvalue <= 0.001){
    star <- "***"
  }
  if(pvalue <= 0.0001){
    star <- "****"
  }
  
  height <- max(samples[, cell_type])
  span <- max(samples[, cell_type]) - min(samples[, cell_type])
  
  plot <- plot +
    annotate("segment", x = 0.9, xend = 2.1, y = (height + 0.03*span), yend = (height + 0.03*span), color = "black") +
    annotate("text", x = 1.5, y = (height + 0.06*span), label = star, size = 5)
  plot
  
  ggsave(paste0("output/Figures/CelltypeProportions/CellTypes_", cell_type, "_CaseCtrl.pdf"), plot, units = "cm", width = 6, height = 9, dpi = 400)
  ggsave(paste0("output/Figures/CelltypeProportions/CellTypes_", cell_type, "_CaseCtrl.svg"), plot, units = "cm", width = 6, height = 9, dpi = 400)
  save(plot, height, span, file = paste0("output/Figures/CelltypeProportions/CellTypes_", cell_type, "_CaseCtrl.RData"))
}


# COMBINED FIGURE
load("output/Figures/CelltypeProportions/CellTypes_Neutro_CaseCtrl.RData")
CellTypesNeutro <- plot
load("output/Figures/CelltypeProportions/CellTypes_Mono_CaseCtrl.RData")
CellTypesMono <- plot
load("output/Figures/CelltypeProportions/CellTypes_CD4T_CaseCtrl.RData")
CellTypesCD4T <- plot
load("output/Figures/CelltypeProportions/CellTypes_NK_CaseCtrl.RData")
CellTypesNK <- plot
load("output/Figures/CelltypeProportions/CellTypes_Epi_CaseCtrl.RData")
CellTypesEpi <- plot
load("output/Figures/CelltypeProportions/CellTypes_B_CaseCtrl.RData")
CellTypesB <- plot
load("output/Figures/CelltypeProportions/CellTypes_Fib_CaseCtrl.RData")
CellTypesFib <- plot


CellTypesNeutro <- CellTypesNeutro + labs(title ="A", subtitle = "Neutrophils") + theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5), plot.title = element_text(hjust = 0))
CellTypesMono <- CellTypesMono + labs(title ="B", subtitle = "Monocytes") + theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5), plot.title = element_text(hjust = 0))
CellTypesCD4T <- CellTypesCD4T + labs(title ="C", subtitle = "CD4T cells") + theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5), plot.title = element_text(hjust = 0))
CellTypesNK <- CellTypesNK + labs(title ="D", subtitle = "NK cells") + theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5), plot.title = element_text(hjust = 0))
CellTypesEpi <- CellTypesEpi + labs(title ="E", subtitle = "Epithelial cells") + theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5), plot.title = element_text(hjust = 0))
CellTypesB <- CellTypesB + labs(title ="F", subtitle = "B cells") + theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5), plot.title = element_text(hjust = 0))
CellTypesFib <- CellTypesFib + labs(title ="G", subtitle = "Fibroblasts") + theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5), plot.title = element_text(hjust = 0))

empty <- ggplot() + theme_void() #+ th + th_transparent
empty2 <- ggplot() + theme_void() #+ th + th_transparent

design <- "
  ABCD
  EEEE
  FGHI
  "

layout <- (CellTypesNeutro + CellTypesMono + CellTypesCD4T + CellTypesNK + empty + CellTypesEpi + CellTypesB + CellTypesFib + empty2) +
  plot_layout(design = design, heights = c(1, 0.1, 1))

ggsave(paste0("output/Figures/CelltypeProportions/Cell_type_prop_categ_COMB.svg"), layout, units = "cm", width = 35, height = 25, dpi = 400)
ggsave(paste0("output/Figures/CelltypeProportions/Cell_type_prop_categ_COMB.pdf"), layout, units = "cm", width = 35, height = 25, dpi = 400)       


### COMBINE PLOTS ###

# 2) VARIATION

celltype_prop_hepi_plot <- celltype_prop_hepi_plot + labs(title ="A") + xlab("Samples") + theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5), plot.title = element_text(hjust = 0))
Hepi_violin_plot <- Hepi_violin_plot + labs(title ="B") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + theme(plot.subtitle = element_text(size = 14, face = "bold", hjust = 0.5), plot.title = element_text(hjust = 0))

design <- "AB"

layout <- (celltype_prop_hepi_plot + Hepi_violin_plot) + plot_layout(design = design)

ggsave(paste0("output/Figures/CelltypeProportions/Variation_COMB.svg"), layout, units = "cm", width = 26, height = 13, dpi = 400)
ggsave(paste0("output/Figures/CelltypeProportions/Variation_COMB.pdf"), layout, units = "cm", width = 26, height = 13, dpi = 400)       


