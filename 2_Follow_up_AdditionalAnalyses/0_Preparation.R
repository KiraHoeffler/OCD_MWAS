
################################################################################
# SETUP
###############################################################################

# PATH TO WORKING DIRECTORY:
dir_gen <- "S://Project/WP-epigenetics/10_CaseCtrl_EWAS/"

################################################################################
# SET UP FOLDER STRUCTURE FOR THE output FILES
################################################################################

# SET WORKING DIRECTORY
setwd(dir_gen)

# output FOLDER STRUCTURE
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/RData")) dir.create("output/RData")
if (!dir.exists("output/Figures")) dir.create("output/Figures")
if (!dir.exists("output/Tables")) dir.create("output/Tables")

#LOAD PACKAGE
library(ggplot2)

# DEFINE THEME
th <- theme(panel.background = element_rect(fill = "white"),
            plot.title = element_text(size = 14, color = "black", face = "bold"),
            axis.line = element_line(colour = "black", linewidth = 0.5),
            axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
            axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
            axis.title.x = element_text(vjust = -0.5, size = 12, colour = "black"),
            axis.title.y = element_text(vjust = -0.5, size = 12, colour = "black", margin = margin(r = 10)),
            axis.title = element_text(face = "bold"),
            panel.border = element_blank(),
            legend.text = element_text(size = 12),
            legend.title = element_text(size =12),
            legend.key = element_blank()) 

th_transparent <- theme(panel.background = element_rect(fill = "transparent"),
                        plot.background = element_rect(fill = "transparent", color = NA),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.background = element_rect(fill = "transparent", color = NA),
                        legend.box.background = element_rect(fill = "transparent"))

save(th, file = "output/RData/theme.RData")
save(th_transparent, file = "output/RData/theme_transparent.Rdata")