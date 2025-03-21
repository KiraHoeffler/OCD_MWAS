
# LOAD LIBRARY
suppressMessages({
  library(readxl)
})


### DECIDE WHICH SCRIPT TO RUN DEPENDING ON THE NR OF MALES AND FEMALES AFTER QC
samplesheet <- as.data.frame(read_xlsx("output/Tables/samplesheet_after_QC.xlsx"))

if (multiple_samples == "yes"){
  samplesheet <- samplesheet[which(samplesheet$Baseline_sample == "Y"),]
}

if (case_ctrl_analysis == "yes"){
  print("Running sex-adjusted case-control analysis.")
  
  source("script/14_CaseControl_Analysis_Females_Males_limma_saliva.R")
  
  # CREATE A REPORT
  print("Create report")
  source("2_Paths.R")
  rmarkdown::render("script/Report_Part4_CaseControl_sex-adjusted.Rmd", 
                    output_file = paste0(dir_gen,"Report_Part5_CaseControl_Analysis_sex_adjusted.html"),
                    params = list(dir_gen = dir_gen)
  )
}