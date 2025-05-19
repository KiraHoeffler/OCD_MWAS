
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
samples_F <- as.data.frame(read_xlsx("output_F/Tables/Final_samplesheet_CaseCtrl.xlsx"))
samples_M <- as.data.frame(read_xlsx("output_M/Tables/Final_samplesheet_CaseCtrl.xlsx"))
samples <- rbind(samples_F, samples_M)
samples <- samples[order(samples$Basename), ]

################################################################################
# SAMPLESHEET
################################################################################

table(samples$Case_Control)
#control case 
#384  414

### OCD ###
samples_OCD <- samples[which(samples$Case_Control == "Case"), ]

#Sex 
table(samples_OCD$Sex)
297/(297+117) #0.7173913

# Age
mean(samples_OCD$Age) #31.45821
sd(samples_OCD$Age) #10.19994

# Ancestry
table(samples_OCD$Ancestry)
table(is.na(samples_OCD$Ancestry))

#Medication
table(samples_OCD$Psychoactive_medicine)
table(is.na(samples_OCD$Psychoactive_medicine))

table(samples_OCD$Antidepressants)
table(samples_OCD$Stimulants)
table(samples_OCD$Antipychotics)
table(samples_OCD$Mood_stablizers)
table(samples_OCD$Sedatives)
table(samples_OCD$Anxiolytics)

#Comorbidities
table(samples_OCD$Comorbidity)
table(samples_OCD$combined_current_depression)
table(samples_OCD$bipolar1_current)
table(samples_OCD$panic_current)
table(samples_OCD$agora_current)
table(samples_OCD$combined_social_phobia_current)
table(samples_OCD$PTSD_current)
table(samples_OCD$anorexia_current_3mo)
table(samples_OCD$bulimia_current_3mo)
table(samples_OCD$GAD_current_6mo)



### CTRL ###
samples_CTRL <- samples[which(samples$Case_Control == "Control"), ]

#Sex 
table(samples_CTRL$Sex)
275/(275+109) #0.7161458

# Age
mean(samples_CTRL$Age) #30.23997
sd(samples_CTRL$Age) #6.812872

# Ancestry
table(samples_CTRL$Ancestry)
table(is.na(samples_CTRL$Ancestry))


samplesheet_old <- as.data.frame(read_xlsx("S:/Project/WP-epigenetics/03_Phenotype_sheet/Archive/samplesheet_v5.xlsx"))
table(is.na(samplesheet_old$anc_classif))

missing_ancestry <- samplesheet_old$Basename[which(is.na(samplesheet_old$anc_classif))]

table(samples_OCD$Basename %in% missing_ancestry)
#FALSE  TRUE 
#397    17 

table(samples_CTRL$Basename %in% missing_ancestry)
#FALSE  TRUE 
#136   248 




# STATISTICAL TESTS

shapiro.test(samples$Age)
#data:  samples$Age
#W = 0.92571, p-value < 2.2e-16 -> not normally distributed

wilcox.test(Age ~ Case_Control, data = samples)
#data:  Age by Case_Control
#W = 78753, p-value = 0.8213
#alternative hypothesis: true location shift is not equal to 0

chisq.test(table(samples$Case_Control, samples$Sex))
#data:  table(samples$Case_Control, samples$Sex)
#X-squared = 5.6569e-30, df = 1, p-value = 1


# Females
 
shapiro.test(samples_F$Age)
#data:  samples_F$Age
#W = 0.9172, p-value < 2.2e-16

wilcox.test(Age ~ Case_Control, data = samples_F)
# data:  Age by Case_Control
# W = 39995, p-value = 0.6698
# alternative hypothesis: true location shift is not equal to 0


# Males

shapiro.test(samples_M$Age)
#data:  samples_M$Age
#W = 0.94028, p-value = 5.472e-08

wilcox.test(Age ~ Case_Control, data = samples_M)
#data:  Age by Case_Control
#W = 6452, p-value = 0.8786
#alternative hypothesis: true location shift is not equal to 0

