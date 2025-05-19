


#' Adapt combp results
#'
#' adapts output of combp function
#' IMPORTANT: THIS FUNCTION WAS WRITTEN FOR EPICv2 - might need some adaptation for other arrays
#' 
#' @param combp_result direct output of combp() function, with columns p.adjust, n, chr, start, end
#' @param DMPs_anno merged DMP results (summary statistics) with annotation, with columns P.Value, CHR, MAPINFO, cgID,
#' @param annotation updated annotation with columns: Gene, RegulRegion, Relation_to_UCSC_CpG_Island
#' @param nr_signif_CpGs DMR should have at least X CpGs (insert a number for X, e.g. 3)
#' @param nr_signif_CpGs DMR should have at least X nominally significant CpGs (insert a number for X, e.g. 3)
#' @param signif_proportion proportion of nominally significant CpGs required in DMR 
#'
#' @return adapted combp df
#'  
adapt_combp <- function(combp_result, DMPs_anno, annotation, nr_CpGs, nr_signif_CpGs, signif_proportion){
  
  print(paste0("significant DMRs before filtering step: ",nrow(combp_result)))
  
  # only select significant DMRs and those with at least 3 CpGs
  combp_sign <- combp_result[which(as.numeric(combp_result$sidak) <= 0.05 & as.numeric(combp_result$nprobe) >= nr_CpGs), ]
  
  # nominal significant DMPs
  DMPs_anno_sign05 <- DMPs_anno[which(DMPs_anno$P.Value < 0.05), ]
  
  #annotate DMP_overlap
  combp_sign$GeneString <- NA
  combp_sign$RegulRegion <- NA
  combp_sign$RelationCpGIsland <- NA
  combp_sign$Nr_CpGs <- NA
  combp_sign$Nr_signif_CpGs_05 <- NA
  
  
  for (i in c(1:nrow(combp_sign))){
    chr <- combp_sign$chr[i]
    start <- combp_sign$start[i]
    end <- combp_sign$end[i]
    
    Genes <- unique(annotation$Gene[which(annotation$CHR == chr & annotation$MAPINFO >= start &  annotation$MAPINFO <= end)])
    Genes <- unique(unlist(strsplit(Genes, ";")))
    GeneString <- paste(Genes, collapse = ";")
    combp_sign$GeneString[i] <- GeneString
    
    RegulRegion <- unique(annotation$RegulRegion[which(annotation$CHR == chr & annotation$MAPINFO >= start &  annotation$MAPINFO <= end)])
    RegulRegion <- unique(unlist(strsplit(RegulRegion, ";")))
    RegulRegionString <- paste(RegulRegion, collapse = ";")
    combp_sign$RegulRegion[i] <- RegulRegionString
    
    RelationCpGIsland <- unique(annotation$Relation_to_UCSC_CpG_Island[which(annotation$CHR == chr & annotation$MAPINFO >= start &  annotation$MAPINFO <= end)])
    RelationCpGIsland <- unique(unlist(strsplit(RelationCpGIsland, ";")))
    RelationCpGIslandString <- paste(RelationCpGIsland, collapse = ";")
    combp_sign$RelationCpGIsland[i] <- RelationCpGIslandString
    
    DMR_cpgs <- unique(annotation$Name[which(annotation$CHR == chr & annotation$MAPINFO >= start & annotation$MAPINFO <= end)])
    combp_sign$Nr_CpGs[i] <- length(DMR_cpgs)
    nr_sign05 <- length(intersect(DMR_cpgs, DMPs_anno_sign05$cgID))
    combp_sign$Nr_signif_CpGs_05[i] <- nr_sign05
  }
  
  combp_sign <- combp_sign[which(combp_sign$Nr_signif_CpGs_05 >= nr_signif_CpGs), ]
  
  combp_sign$Perc_sign_CpGs <- combp_sign$Nr_signif_CpGs_05 / (combp_sign$Nr_CpGs)
  combp_sign <- combp_sign[which(combp_sign$Perc_sign_CpGs >= signif_proportion), ]
  
  print(paste0("significant DMRs after filtering: ", nrow(combp_sign)))
  
  # oder: adj. p value
  combp_sign <- combp_sign[order(combp_sign$sidak), ]
  
  # add UNANNOTATED
  length_unannotated <- length(combp_sign$GeneString[which(combp_sign$GeneString == "")])
  
  count <- 1
  for (i in c(1:nrow(combp_sign))){
    if (combp_sign$GeneString[i] == ""){
      New_Gene <- paste0("UNANNOTATED", count)
      combp_sign$GeneString[i] <- New_Gene
      count <- count + 1
    }
  }
  
  return(combp_sign)
}







#' ADD LOCATION RELATION TO GENE & CPG ISLAND
#'
#' @param combp_sign output of adapt_combp
#' @param annotation  annotation as downloaded from Illumina (EPICv2)
#'
#' @return adapted combp_sign
#'
add_location <- function(combp_sign, annotation){
  
  #annotate DMP_overlap
  combp_sign$Gene_relation <- NA
  combp_sign$CpG_island_relation <- NA
  
  for (i in c(1:nrow(combp_sign))){
    chr <- combp_sign$chr[i]
    start <- combp_sign$start[i]
    end <- combp_sign$end[i]
    
    Gene_anno <- unique(annotation$UCSC_RefGene_Group[which(annotation$CHR == chr & annotation$MAPINFO >= start & annotation$MAPINFO <= end)])
    Gene_relation <- unique(c(unlist(strsplit(Gene_anno, ";"))))
    Gene_relation_String <- paste(sort(Gene_relation), collapse = ";")
    combp_sign$Gene_relation[i] <- Gene_relation_String
    
    Gene_Island <- unique(annotation$Relation_to_UCSC_CpG_Island[which(annotation$CHR == chr & annotation$MAPINFO >= start & annotation$MAPINFO <= end)])
    GeneIsland_relation <- unique(c(unlist(strsplit(Gene_Island, ";"))))
    GeneIsland_relation_String <- paste(sort(GeneIsland_relation), collapse = ";")
    combp_sign$CpG_island_relation[i] <- GeneIsland_relation_String
    
  }
  
  return(combp_sign)
}


###################



#' Extract CpGs from DMRs
#'
#' extract all CpG names from the DMRs
#' 
#' IMPORTANT: written for EPICv2 - might need some adaptation to use for the other arrays
#' 
#' @param combp_sign output of the adapt_combp function
#' @param anno annotation file containing only the columns c("Name", "CHR", "MAPINFO")
#'
#' @return vector with CpGs extracted from DMRs
#' 
extract_CpGs_from_DMRs <- function(combp_sign, anno){
  all_DMR_cpgs <- c()
  
  for (i in c(1:nrow(combp_sign))){
    chr <- combp_sign$chr[i]
    start <- combp_sign$start[i]
    end <- combp_sign$end[i]
    DMR_cpgs <- unique(anno$Name[which(anno$CHR == chr & anno$MAPINFO >= start & anno$MAPINFO <= end)])
    all_DMR_cpgs <- c(all_DMR_cpgs, DMR_cpgs)
  }
  
  return(all_DMR_cpgs)
}



#####################


#' Minimum p in DMR
#' 
#' IMPORTANT: written for EPICv2 - might need some adaptation to use for the other arrays
#'
#' @param DMPs_anno merged DMP results (summary statistics) with annotation, with columns P.Value, cgID
#' @param all_DMR_cpgs output from extract_CpGs_from_DMRs() function
#'
#' @return lowest p value in DMR (needed for Manhattan plot)
#'
min_p_in_DMR <- function(DMPs_anno, all_DMR_cpgs){
  all_DMR_cpgs_pvalues <- DMPs_anno$P.Value[which(DMPs_anno$cgID %in% all_DMR_cpgs)]
  min_DMR_p <- -log10(min(all_DMR_cpgs_pvalues))
  return(min_DMR_p)
}


###################



#' Manhattan plot
#' 
#' IMPORTANT: written for EPICv2 - might need some adaptation to use for the other arrays
#'
#' @param DMPs_anno merged DMP results (summary statistics) with annotation, with columns P.Value, CHR, MAPINFO, cgID
#' @param sex "F" or "M" or "MF" (both sexes), "FM_auto" (only autosomes)
#' @correction "BH" or "bonferroni"  
#' @param combp_sign output from adapt_combp() function
#' @min_DMR_p lowest p value of any DMP in any DMR
#' @min_sign max sign p value (BH threshold)
#'
#' @return Manhatten plot object
#' 
Manhattan <- function(DMPs_anno, sex, correction, combp_sign){
  
  if (correction == "BH"){
    signif_DMPs <- DMPs_anno[which(DMPs_anno$adj_p_BH <= 0.05),]
  }
  if (correction == "bonferroni"){
    signif_DMPs <- DMPs_anno[which(DMPs_anno$adj_p_BH <= 0.05),]
  }
  
  # EXTRACT CPGS IN DMRs
  all_DMR_cpgs <- extract_CpGs_from_DMRs(combp_sign, anno_upd)
  
  # MINIMUM P VALUE IN DMR
  min_DMR_p <- min_p_in_DMR(DMPs_anno, all_DMR_cpgs)
  
  min_sign <- -log10(max(signif_DMPs_BH$P.Value))
  
  
  # ADD -LOG10 P VALUES, MAPINFO, AND COLOR TO THE RESULTS
  DMPs_anno$pvalue_Manhattan <- -log10(DMPs_anno$P.Value)
  DMPs_anno$color <- NA
  DMPs_anno$CHR <- substr(DMPs_anno$CHR, 4, nchar(DMPs_anno$CHR))
  odd <- c("1","3","5","7","9","11","13","15","17","19","21","X")
  DMPs_anno$color <- ifelse(DMPs_anno$CHR %in% odd, "#d4e8e6", "#a8d1ce")
  DMPs_anno$color[which(DMPs_anno$cgID %in% all_DMR_cpgs)] <- "#3c726d"
  DMPs_anno$color[which(DMPs_anno$cgID %in% signif_DMPs$cgID)] <- "#192f2d"
  DMPs_anno$chr2 <- paste0("chr", DMPs_anno$CHR)
  
  if (sex == "F" | sex == "FM"){
    DMPs_anno$chr2 <- factor(DMPs_anno$chr2, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                                        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                                                        "chr19", "chr20", "chr21", "chr22", "chrX"))
  }
  
  if (sex == "FM_auto"){
    DMPs_anno$chr2 <- factor(DMPs_anno$chr2, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                                        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                                                        "chr19", "chr20", "chr21", "chr22"))
  }
  
  if (sex == "M"){
    DMPs_anno$chr2 <- factor(DMPs_anno$chr2, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                                        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                                                        "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))
  }
  
  #continuous counting of positions
  count_data <- DMPs_anno %>%
    group_by(chr2) %>%
    summarise(max_position = max(MAPINFO)) %>%
    mutate(position_add = lag(cumsum(as.numeric(max_position)), default = 0)) %>%
    select(chr2, position_add)
  
  #add actual position with continous chromosome counting position
  DMPs_anno2 <- DMPs_anno %>%
    inner_join(count_data, by = "chr2") %>%
    mutate(position_cum = MAPINFO + position_add)
  
  #position for x axis chr label
  axis_point <- DMPs_anno2 %>%
    group_by(chr2) %>%
    summarise(center = mean(position_cum))
  
  if (sex == "F" | sex == "FM"){
    axis_point$CHR <- c("1","2","3","4","5","6","7","8","","10","","12","","14","","16","","18","","","","", "X")
  }
  
  if (sex == "FM_auto"){
    axis_point$CHR <- c("1","2","3","4","5","6","7","8","","10","","12","","14","","16","","18","","","","")
  }
  
  
  if (sex == "M"){
    axis_point$CHR <- c("1","2","3","4","5","6","7","8","","10","","12","","14","","16","","18","","","","", "X", "Y")
  }
  
  
  # y limit
  y_limit <- max(DMPs_anno2$pvalue_Manhattan) + 0.5
  
  # DMR coordinates
  DMR_coordinates <- DMPs_anno2$position_cum[which(DMPs_anno2$color == "#3c726d")]
  
  #Bonferroni
  #min_p <- 0.05/nrow(DMPs_anno)
  
  # ACTUAL PLOT
  Manhattan_plot <- ggplot() +
    geom_point(data = subset(DMPs_anno2, color == "#d4e8e6"), aes(x=position_cum, y=pvalue_Manhattan, color = color), shape = 19, size = 0.8) +
    geom_point(data = subset(DMPs_anno2, color == "#a8d1ce"), aes(x=position_cum, y=pvalue_Manhattan, color = color), shape = 19, size = 0.8) +
    geom_point(data = subset(DMPs_anno2, color == "#3c726d"), aes(x=position_cum, y=pvalue_Manhattan, color = color), shape = 19, size = 0.8) +
    geom_point(data = subset(DMPs_anno2, color == "#192f2d"), aes(x=position_cum, y=pvalue_Manhattan, color = color), shape = 19, size = 0.8) +
    scale_size_manual(values = c(1,2)) +
    th + th_transparent +
    theme(axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(size = 10)) +
    labs(x = "genomic position (chromosome)", y = "-log10(p)") +
    scale_x_continuous(expand = c(0,0), label = axis_point$CHR, breaks = axis_point$center) + 
    scale_y_continuous(expand = c(0,0), limits = c(0, y_limit)) +
    scale_color_identity() +
    geom_hline(yintercept = 5, color = "grey", linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = min_sign, color = "grey", linetype = "dashed", alpha = 0.5) +
    lapply(DMR_coordinates, function(x) geom_segment(aes(x = x, xend = x, y = 0, yend = min_DMR_p), linewidth = 0.25, color = "#3c726d"))
  
  return(Manhattan_plot)
}


###################

add_avgDNAm <- function(samplesheet, DNAm_DMR, IndivName, cols){
  
  samplesheet_red <- unique(samplesheet[, c(IndivName, cols)])
  samplesheet_red$AvgDMRmethylation <- NA
  
  for (i in c(1:nrow(samplesheet_red))){
    Indiv_name <- samplesheet_red[i, IndivName]
    samples <- samplesheet$Basename[which(samplesheet[, IndivName] == Indiv_name)]
    if (class(beta_DMR)[1] == "numeric"){
      DNAm_ind <- DNAm_DMR[samples]
    } else {
      DNAm_ind <- DNAm_DMR[, samples]
    }
    samplesheet_red$AvgDMRmethylation[i] <- mean(DNAm_ind)
  }
  
  samplesheet_red <- na.omit(samplesheet_red)
  return(samplesheet_red)
}



#################

DMR_corr_plot <- function(DMR_contin, GeneString, corr){
  
  Corr_plot <- ggplot(DMR_contin, aes(x = Percent_improvement, y = AvgDMRmethylation)) +
    geom_point(color = "#084081", alpha = 0.15) + 
    geom_smooth(color = "#084081", method = "lm", linewidth=1.2, alpha = 0.35) +
    scale_color_identity() +
    ggtitle(GeneString) +
    th + xlab("% clinical improvement") + ylab("Average DNAm of CpGs in DMR") + th_transparent +
    theme(plot.title = element_text(hjust = 0.5)) + 
    annotate("text", x=15, y = min(DMR_contin$AvgDMRmethylation), label = paste0("r = ", corr), color = "black")
  return(Corr_plot)
}

DMR_corr_plot2 <- function(DMR_contin, GeneString, corr){
  
  Corr_plot <- ggplot(DMR_contin, aes(x = Percent_improvement, y = AvgDMRmethylation)) +
    geom_point(color = "#084081", alpha = 0.15) + 
    geom_smooth(color = "#084081", method = "lm", linewidth=1.2, alpha = 0.35) +
    scale_color_identity() +
    ggtitle(GeneString) +
    th + xlab("% clinical improvement") + ylab("Average DNAm of CpGs in DMR") + th_transparent +
    theme(plot.title = element_text(hjust = 0.5)) + 
    annotate("text", x=20, y = min(DMR_contin$AvgDMRmethylation), label = paste0("r = ", corr), color = "black")
  return(Corr_plot)
}


##################

DMR_categ_plot <- function(DMR_categ, GeneString){
  Categ_plot <- ggplot(DMR_categ, aes(x = Response_status, y = AvgDMRmethylation, fill = Response_status)) +
    geom_violin(scale = "width", alpha = 0.3) + 
    geom_boxplot(width = 0.1, fill = "white") +
    geom_point(data = df_summary, aes(y = Mean, x = Response_status), color = "#14747b", size = 1) +
    th + th_transparent +
    ggtitle(GeneString) +
    labs(x = NULL, y = "Average DNAm of CpGs in DMR") + 
    scale_fill_manual(values = c("resp" = "#084081", "non-resp" = "#891a1a", "ctrl" = "#7bccc4"))+
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
  return(Categ_plot)
}


###############


#' Plot a region with combp.plot
#' -> adapted original function, just changed the look a bit
#'
#' Calculate statistics for a set of genomic regions.
#'
#' @param region.chr/start/end Genomic region to plot.
#' @param estimate Vector of EWAS effect estimates
#' (corresponds to rows of \code{methylation}).
#' @param se Vector of standard errors of the coefficients.
#' @param methylation Methylation matrix (rows=features, columns=samples).
#' @param chr Feature chromosome (corresponds to rows of \code{methylation}).
#' @param pos Feature chromosome position.
#' @param expand Proportion to expand to include sites beyond the beginning and end of the region (default: 1, i.e. plot 3x the length of the region).
#' @param ci Show confidence intervals rather than standard errors (default: TRUE).
#' @param verbose If \code{TRUE} (default), then output status messages.
#' @return \code{\link{ggplot}} object showing the plot. 
#'
#' @export
combp.plot <- function(region.chr, region.start, region.end, estimate, se, chr, pos, expand=1, ci=T, verbose=T) {
  require(ggplot2)
  stopifnot(length(region.chr)==1)
  stopifnot(length(region.start)==1)
  stopifnot(length(region.end)==1)
  stopifnot(is.numeric(region.start))
  stopifnot(is.numeric(region.end))
  stopifnot(region.start < region.end)
  stopifnot(is.vector(estimate))
  stopifnot(is.vector(se))
  stopifnot(is.vector(chr))
  stopifnot(is.vector(pos))
  stopifnot(length(estimate)==length(se))
  stopifnot(length(estimate)==length(chr))
  stopifnot(length(estimate)==length(pos))
  stopifnot(region.chr %in% chr)
  region.len <- region.end-region.start
  plot.start <- max(0,region.start - region.len*expand)
  plot.end <- region.end + region.len*expand
  idx <- which(chr==region.chr & pos >= plot.start & pos <= plot.end)
  if (length(idx) == 0)
    stop("The region does not contain any CpG sites")
  dat <- data.frame(estimate,se,pos)[idx,]
  dat$ymax <- dat$estimate + dat$se*ifelse(ci,1.96,1)
  dat$ymin <- dat$estimate - dat$se*ifelse(ci,1.96,1)
  dat <- dat[order(dat$pos),]
  gap <- max(1, min(tail(dat$pos,-1)-head(dat$pos,-1), na.rm=T)/2)
  (ggplot(dat, aes(x=pos, y=estimate, ymin=ymin, ymax=ymax)) +
      geom_rect(
        aes(xmin=region.start-gap, xmax=region.end+gap, ymin=-Inf, ymax=Inf),
        fill="#BBBBBB", alpha=0.1) +
      ggtitle(GeneString) +
      geom_hline(yintercept=0, linetype="dashed", color="black") +
      geom_pointrange() +
      labs(
        y=paste0("estimate (", ifelse(ci, "95% CI", "SE"), ")"),
        x=paste0(region.chr, " position")) +
      th + th_transparent +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))    
}
