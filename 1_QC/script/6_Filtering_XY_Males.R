
################################################################################
# SETUP
################################################################################

print("6 - FILTERING X AND Y CHR IN MALES")
print("Setting up and loading libraries")

# LOAD PATHS
source("2_Paths.R")


### OUTLIER SELECTION (ADD OUTLIER WHEN THERE IS STILL AN OUTLIER AFTER RUNNING THE SCRIPT FOR THE FIRST TIME, SEE GUIDELINES) ###

outlier <- c() #identify outlier at the end of the script

########################

# SET DETECTION P VALUE THRESHOLD
thres=1E-5 #detection p-value threshold

if (array.type == "EPICv1"){
  annofile_path <- file.path(resources,"EPICv1_manifest.csv") #downloaded Illumina annotation file
  probes_to_excl_path <- file.path(resources,"EPICv1_RemoveProbes.Robj")
  SNP_probes_path <- file.path(resources,"EPICv1_SNP_probes_0bp_MAF05.RData")
}
if (array.type == "EPICv2"){
  annofile_path <- file.path(resources,"EPICv2_manifest.csv")
  probes_to_excl_path <- file.path(resources,"EPICv2_RemoveProbes.RData")
  SNP_probes_path <- file.path(resources,"EPICv2_SNP_probes_0bp_MAF05_IlmnID.RData")
}
if(array.type == "450K"){
  annofile_path <- file.path(resources,"450k_manifest.csv")
  probes_to_excl_path <- file.path(resources,"450K_RemoveProbes_ChenProbeIDs.Robj")
  SNP_probes_path <- file.path(resources,"450K_SNP_probes_0bp_MAF05.RData")
}

callrate_thresh <- 0.95 # call rate threshold

################################################################################
# LOAD LIBERIES
################################################################################

suppressMessages({
  #library(DBI, lib.loc = "Z:/Bioconductorpackages/319")
  #library(tidyselect, lib.loc = "Z:/Bioconductorpackages/319")
  
  # LOAD LIBRARIES
  library(minfi)
  library(readxl)
  library(ENmix)
  library(tidyr)
  library(ggplot2)
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

# SET WORKING DIRECTORY
setwd(dir_gen)

# LOAD INTENSITIES & DETECTION P VALUES
load('output/RData/intensities.RData')
load('output/RData/detectionPvalue.RData')

# IMPORT ANNOTATION FILE
annotation <- read.csv(annofile_path, skip = 7)

# LOW BEAD COUNT
load("output/RData/Low_beadcount.RData")

# IMPORT SAMPLESHEET
samplesheet <- as.data.frame(read_excel("output/Tables/Samplesheet.xlsx"))

# PROBES TO EXCLUDE
load(probes_to_excl_path)


################################################################################
# GET CPGS ON XY & SAMPLE NAMES OF MALES
################################################################################

# EXTRACT CPGs ON SEX CHROMOSOMES FROM ANNOTATION FILE
if (array.type == "EPICv2"){
  cgs <- annotation[substr(annotation$Name, start = 1, stop = 2) == "cg" & (annotation$CHR == "chrX" | annotation$CHR == "chrY"),]
}

if (array.type == "EPICv1" | array.type =="450K"){
  cgs <- annotation[substr(annotation$Name, start = 1, stop = 2) == "cg" & (annotation$CHR == "X" | annotation$CHR == "Y"),]
}


auto <- as.matrix(cgs$IlmnID)
rm(annotation)

# SAMPLE NAMES
samples <- intersect(colnames(TypeI.Red.M), samplesheet$Basename[which(samplesheet$Sex == "M")])


################################################################################
# APPLY DETECTION P VALUE THRESHOLD
################################################################################

print("Step 1/5: apply detection p value threshold")

# SET VALUES ABOVE DETECTION P VALUE THRESHOLD TO NA IN INTENSITIES DFS
d=dp[rownames(TypeII.Green),colnames(TypeII.Green)]
TypeII.Green.d = ifelse(d<thres,TypeII.Green,NA)
TypeII.Red.d = ifelse(d<thres,TypeII.Red,NA)
d=dp[rownames(TypeI.Green.M),colnames(TypeI.Green.M)]
TypeI.Green.M.d = ifelse(d<thres,TypeI.Green.M,NA)
TypeI.Green.U.d = ifelse(d<thres,TypeI.Green.U,NA)
d=dp[rownames(TypeI.Red.M),colnames(TypeI.Red.M)]
TypeI.Red.M.d = ifelse(d<thres,TypeI.Red.M,NA)
TypeI.Red.U.d = ifelse(d<thres,TypeI.Red.U,NA)
rm(dp,d)

################################################################################
# FILTER INTENSITIES ON SAMPLES & X CpGs
################################################################################

print("Step 2/5: filter intensities on XY chr CpGs and males")

markers=as.matrix(intersect(rownames(TypeII.Green.d), auto))
TypeII.Green = TypeII.Green.d[markers,samples]
TypeII.Red = TypeII.Red.d[markers,samples]
markers=intersect(rownames(TypeI.Green.M.d), auto)
TypeI.Green.M = TypeI.Green.M.d[markers,samples]
TypeI.Green.U = TypeI.Green.U.d[markers,samples]
markers=intersect(rownames(TypeI.Red.M.d), auto)
TypeI.Red.M = TypeI.Red.M.d[markers,samples]
TypeI.Red.U = TypeI.Red.U.d[markers,samples]



################################################################################
# SAMPLES AND PROBES TO EXCLUDE
################################################################################

print("Step 3/5: Exclude samples and probes that do not pass the QC.")

# SEX  MISMATCH
samples_to_excl <- c(samplesheet$Basename[which(samplesheet$Sex_mismatch == "Sex_mismatch" | samplesheet$Sex_outlier == "Sex_outlier") ])
samples_to_excl <- unique(c(samples_to_excl, outlier))

samplesheet$to_exclude1e16_XY <- NA
samplesheet$to_exclude1e16_XY[which(samplesheet$Basename %in% samples_to_excl)] <- "exclude"

# HANDLE ADDITIONAL OUTLIERS
if (length(additional_outliers) > 0){
  samples_to_excl <- unique(c(samples_to_excl, additional_outliers))
}

# COMBINE
probes_to_excl <- unique(c(remove_probes, lowbeadcount_df$CpG)) 

################################################################################
# FILTERED INDIVIDUALS / MARKERS
################################################################################

# FILTERED SAMPLES/INDIVIDUALS
auto_red <- setdiff(auto, probes_to_excl) #to keep
samples_red <- setdiff(samples, samples_to_excl)#to keep

# FILTER ON PROBES AND INDIVIDUALS
markers=intersect(rownames(TypeII.Green.d), auto_red)
TypeII.Green = TypeII.Green.d[markers,samples_red]
TypeII.Red = TypeII.Red.d[markers,samples_red]
markers=intersect(rownames(TypeI.Green.M.d), auto_red)
TypeI.Green.M = TypeI.Green.M.d[markers,samples_red]
TypeI.Green.U = TypeI.Green.U.d[markers,samples_red]
markers=intersect(rownames(TypeI.Red.M.d), auto_red)
TypeI.Red.M = TypeI.Red.M.d[markers,samples_red]
TypeI.Red.U = TypeI.Red.U.d[markers,samples_red]


################################################################################
# FILTERED BETA VALUES
################################################################################

# CALCULATE BETAS FOR DIFFERENT PROBE TYPES & COMBINE
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))


################################################################################
# CALCULATE SAMPLE & MARKER (PROBE) CALL RATES ON FILERED BETA VALUES
################################################################################

print("Step 4/5: exclude samples and probes with low call rates.")

sample.call=colSums(!is.na(beta))/nrow(beta)
marker.call=rowSums(!is.na(beta))/ncol(beta)

# IDENTIFY CPGs WITH A LOW MARKER CALL RATE & EXPORT
CpGs_low_call_rate <- names(marker.call[marker.call < callrate_thresh])
save(CpGs_low_call_rate, file = "output/RData/Low_marker_call_rate_1e16_XY_males.RData")

# IDENTIFY SAMPLES WITH A LOW SAMPLE CALL RATE & EXPORT
low_sample_call_rate <- sample.call[sample.call < callrate_thresh]
samples_low_call_rate <- names(sample.call[sample.call < callrate_thresh])
save(samples_low_call_rate, file = "output/RData/Low_sample_call_rate_1e16_XY_males.RData")
low_sample_call_df <- data.frame (sample = samples_low_call_rate, call_rate_1e16 = low_sample_call_rate)
save(sample.call, marker.call, file='output/RData/callRates_1e16_XY_males.RData')

# MERGE LOW SAMPLE CALL RATE WITH SAMPLESHEET
samplesheet$Low_sample_call_rate[which(samplesheet$Basename %in% low_sample_call_df$sample)] <- "Y"

samplesheet$to_exclude1e16_XY[which(samplesheet$Basename %in% low_sample_call_df$sample)] <- "exclude"
samplesheet$to_exclude1e16_XY[which(samplesheet$Basename %in% outlier)] <- "exclude"

# ADD TO PROBES/SAMPLES TO EXCLUDE
probes_to_excl <- unique(c(probes_to_excl, CpGs_low_call_rate))
all_exclude <- samplesheet$Basename[which(samplesheet$to_exclude1e16 == "exclude")]
samples_to_excl <- unique(c(samples_to_excl, samples_low_call_rate, all_exclude))
print(paste0("Nr of samples to excl: ", length(samples_to_excl)))

################################################################################
# FILTERED INDIVIDUALS / MARKERS
################################################################################

print("Step 5/5: save filtered intensities and visualise filtered DNA methylation distributions.")

# FILTERED SAMPLES/INDIVIDUALS TO KEEP
auto_red <- setdiff(auto, probes_to_excl)
samples_red <- setdiff(samples, samples_to_excl) 
print(paste0("Nr of samples to keep: ", length(samples_red)))

# FILTER ON PROBES AND INDIVIDUALS
markers=intersect(rownames(TypeII.Green.d), auto_red)
TypeII.Green = TypeII.Green.d[markers,samples_red]
TypeII.Red = TypeII.Red.d[markers,samples_red]
markers=intersect(rownames(TypeI.Green.M.d), auto_red)
TypeI.Green.M = TypeI.Green.M.d[markers,samples_red]
TypeI.Green.U = TypeI.Green.U.d[markers,samples_red]
markers=intersect(rownames(TypeI.Red.M.d), auto_red)
TypeI.Red.M = TypeI.Red.M.d[markers,samples_red]
TypeI.Red.U = TypeI.Red.U.d[markers,samples_red]
save(TypeII.Red ,TypeII.Green ,TypeI.Red.M ,TypeI.Red.U ,TypeI.Green.M ,TypeI.Green.U , 
     file="output/RData/intensities_filtered_1e16_XY_males.RData",compress=F)

# CALCULATE BETAS FOR DIFFERENT PROBE TYPES & COMBINE
TypeII.betas = TypeII.Green/(TypeII.Red+TypeII.Green+100)
TypeI.Green.betas = TypeI.Green.M/(TypeI.Green.M+TypeI.Green.U+100)
TypeI.Red.betas = TypeI.Red.M/(TypeI.Red.M+TypeI.Red.U+100)
beta = as.matrix(rbind(TypeII.betas,TypeI.Green.betas,TypeI.Red.betas))

# REPLACE ALL ZEROS IN BETA VALUE WITH LOWEST VALUE
replacement_zero <- min(beta[which(beta > 0 & !is.na(beta))]) 
beta_temp <- beta
beta_temp[beta_temp == 0] <- replacement_zero

# CALCULATE M VALUES 
Mvalues <- B2M(beta_temp)
rm(beta_temp)


# EXPORT
#save(beta, file='output/RData/beta_filtered_1e16.RData',compress=F)
#save(Mvalues, file='output/RData/Mvalues_filtered_1e16.RData',compress=F)

# EXPORT SAMPLESHEET
write_xlsx(samplesheet, path = "output/Tables/Samplesheet.xlsx")

################################################################################
# VISUALISATION OF FILTERED DATA
################################################################################

# VISUALISE FILTERED BETA DISTRIBUTION
pdf("output/Figures/Distributions/Beta_distribution_filtered_1e16_XY_males.pdf", width = 7, height = 5)
densityPlot(beta, main = "Filtered beta values")
dev.off()

# VISUALISE FILTERED M VALUE DISTRIBUTION
pdf("output/Figures/Distributions/Mvalue_distribution_filtered_1e16_XY_males.pdf", width = 7, height = 5)
densityPlot(Mvalues, main = "Filtered M values")
dev.off()



################################################################################
# REMAINING OUTLIER
################################################################################

dens_df <- as.data.frame(matrix(data = NA, nrow = 50, ncol = ncol(beta)))
colnames(dens_df) <- colnames(beta)

for (i in c(1:ncol(beta))){
  dens <- density(beta[,i], na.rm = TRUE, n=50)
  dens_y <- dens$y
  dens_df[, i] <- dens_y
}

dens_df_long <- gather(dens_df, key = "Sample", value = "Value")
dens_df_long$X <- rep(1:50, times = ncol(beta))

suppressWarnings({
  dens_plot_lowintensity <- ggplot(dens_df_long, aes(x = X, y = Value, group = Sample)) + #color = (Sample %in% samples_low_intensity),
    geom_line(size = 0.3) +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "grey")) +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = seq(min(dens_df_long$X), max(dens_df_long$X), by = 2), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, max(dens_df_long$Value), by = 0.2), expand = c(0,0))
  dens_plot_lowintensity
})

ggsave("output/Figures/Distributions/Filtered_beta_outlier_check_XY_males.pdf", plot = dens_plot_lowintensity, device = "pdf", width=10, height=6, dpi=400)

### OUTLIERS ###

# if there is an outlier in the filtered beta distribution plot, select an X value,
# and Y values that separate the outlier from the other samples. Run the code below to identify the outlier:

# RUN THIS CODE (REMOVE THE # AT THE START OF THE NEXT TWO LINES)
#outlier <- dens_df_long$Sample[which(dens_df_long$X == x_value & dens_df_long$Value <= y_upper_limit & dens_df_long$Value >= y_lower_limit)]
#print(outlier)


rm(list = ls());gc()
