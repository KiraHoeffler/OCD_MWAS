
################################################################################
# SET-UP
################################################################################

# SET WORKING DIRECTORY
setwd("S:/Project/WP-epigenetics/10_CaseCtrl_EWAS/")

# LOAD PACKAGES
library(readxl)
library(writexl)

################################################################################
# SAMPLESHEET
################################################################################

# IMPORT SAMPLESHEET
samplesheet <- as.data.frame(read_xlsx("S:/Project/WP-epigenetics/04_Pipeline_PhaseI/output/Tables/Sex_adjusted/Final_samplesheet_CaseCtrl.xlsx"))

################################################################################
# SAMPLESHEET
################################################################################

table(samplesheet$Case_Control)
#control non-responder     responder 
#384            75           329

### OCD ###
samplesheet_OCD <- samplesheet[which(samplesheet$Case_Control == "Case"), ]

#Sex 
table(samplesheet_OCD$Sex)
297/(297+117) #0.7173913

# Age
mean(samplesheet_OCD$Age) #31.45821
sd(samplesheet_OCD$Age) #10.19994

# Ancestry
table(samplesheet_OCD$Ancestry)
table(is.na(samplesheet_OCD$Ancestry))

#Medication
table(samplesheet_OCD$Psychoactive_medicine)
table(is.na(samplesheet_OCD$Psychoactive_medicine))

table(samplesheet_OCD$Antidepressants)
table(samplesheet_OCD$Stimulants)
table(samplesheet_OCD$Antipychotics)
table(samplesheet_OCD$Mood_stablizers)
table(samplesheet_OCD$Sedatives)
table(samplesheet_OCD$Anxiolytics)

#Comorbidities
table(samplesheet_OCD$Comorbidity)
table(samplesheet_OCD$combined_current_depression)
table(samplesheet_OCD$bipolar1_current)
table(samplesheet_OCD$panic_current)
table(samplesheet_OCD$agora_current)
table(samplesheet_OCD$combined_social_phobia_current)
table(samplesheet_OCD$PTSD_current)
table(samplesheet_OCD$anorexia_current_3mo)
table(samplesheet_OCD$bulimia_current_3mo)
table(samplesheet_OCD$GAD_current_6mo)



### CTRL ###
samplesheet_CTRL <- samplesheet[which(samplesheet$Case_Control == "Control"), ]

#Sex 
table(samplesheet_CTRL$Sex)
275/(275+109) #0.7161458

# Age
mean(samplesheet_CTRL$Age) #30.23997
sd(samplesheet_CTRL$Age) #6.812872

# Ancestry
table(samplesheet_CTRL$Ancestry)
table(is.na(samplesheet_CTRL$Ancestry))


samplesheet_old <- as.data.frame(read_xlsx("S:/Project/WP-epigenetics/03_Phenotype_sheet/Archive/samplesheet_v5.xlsx"))
table(is.na(samplesheet_old$anc_classif))

missing_ancestry <- samplesheet_old$Basename[which(is.na(samplesheet_old$anc_classif))]

table(samplesheet_OCD$Basename %in% missing_ancestry)
table(samplesheet_CTRL$Basename %in% missing_ancestry)





# STATISTICAL TESTS

shapiro.test(samplesheet$Age)
#data:  samplesheet$Age
#W = 0.92571, p-value < 2.2e-16

wilcox.test(Age ~ Case_Control, data = samplesheet)
#data:  Age by Case_Control
#W = 78753, p-value = 0.8213
#alternative hypothesis: true location shift is not equal to 0

chisq.test(table(samplesheet$Case_Control, samplesheet$Sex))
#data:  table(samplesheet$Case_Control, samplesheet$Sex)
#X-squared = 5.6569e-30, df = 1, p-value = 1

