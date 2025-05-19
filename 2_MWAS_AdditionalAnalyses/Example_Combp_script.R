
load("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/DMPs/CaseCtrl_Females/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_anno_BACON.RData")
setwd("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/DMPs/CaseCtrl_Females/")


library(ENmix)


combp_input <- data.frame(chr = DMPs_anno$CHR, start = DMPs_anno$MAPINFO, end = DMPs_anno$MAPINFO + 1, 
                          p = DMPs_anno$P.Value, probe = DMPs_anno$cgID)

combp(combp_input, dist.cutoff = 750, seed = 0.001, region_plot = F, mht_plot = F)
