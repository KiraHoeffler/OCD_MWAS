

#####################################################################################################
# Update below (no need to run the script, just adapt and save)
#####################################################################################################

# IS YOUR COHORT INCLUDED IN THE CASE-CONTROL ANALYSIS? ("yes" or "no")
case_ctrl_analysis <- "yes"

# MULTIPLE SAMPLES PER INDIVIDUAL ("yes" or "no")
multiple_samples <- "yes" #then please have a column Baseline_sample with entry "Y"
                          #for all samples that should be included in the analysis

# WHICH ARRAY VERSION DID YOU USE ("450K", "EPICv1", "EPICv2")
array.type <- "EPICv2"

# WHICH TISSUE TYPE ("saliva" or "blood")
tissue_type <- "saliva"

# PATH TO WORKING DIRECTORY (HAS TO CONTAIN THE "script" FOLDER):
dir_gen <- "S:/Project/WP-epigenetics/04_Pipeline_PhaseI/" #working directory

# PATH TO RESOURCES (PROVIDED BY US):
resources = "S:/Project/WP-epigenetics/02_Import/OCD_Pipeline/Resources/"

# PATH TO SAMPLESHEET (details on the samplesheet in our guide)
samplesheet_path <- "S:/Project/WP-epigenetics/04_Pipeline_PhaseI/samplesheet_PhaseI.xlsx"

# PATH TO THE IDAT FILES (RAW DNA METHYLATION DATA)
#dir_idat <- "S:/Project/WP-epigenetics/01_Raw_data/2023_April/2022-086-ILL_GSAUHB_METUHB_N=1672/2022-086-ILL_GSAUHB_METUHB_N=1672/UHB/"?
dir_idat <- "S:/Project/WP-epigenetics/01_Raw_data/2023_April/2022-086-ILL_GSAUHB_METUHB_N=1672/2022-086-ILL_GSAUHB_METUHB_N=1672/UHB/"


# DO YOU HAVE GENOTYPING DATA FOR ALL SAMPLES IN THE SAMPLESHEET ("yes" or "no")
available_genotyping_data <- "no"

if (available_genotyping_data == "yes"){
  # add the path to the quality controlled genotyping data that can be used to calculate PCs
  # please have the data saved in a matrix called genotype_matrix as RData object (samples as columns, SNPs as rows.
  genotyping_path <- "add_complete_path"
}

########## # ADDITIONAL OUTLIERS (OPTIONAL) #################

# if you want to exclude additional samples in the filtering step for whatever reason, please put their Basename here, you can also add a reason, e.g. "genotype_mismatch"
# one reason for each outlier, so the number of outliers and reasons need to be the same

additional_outliers <- c("207036380041_R05C01", "207039770069_R03C01", "207049880010_R02C01", "207057140076_R01C01",
                         "207049880149_R08C01", "207036380163_R01C01", "207036380163_R04C01", "207036380163_R02C01")
reason_exclusion <- c("Genotype_mismatch", "Genotype_mismatch", "Genotype_mismatch", "Genotype_mismatch",
                      "Female_sex_outlier", "Male_sex_outlier", "Male_sex_outlier", "Male_sex_outlier", "Male_sex_outlier")
