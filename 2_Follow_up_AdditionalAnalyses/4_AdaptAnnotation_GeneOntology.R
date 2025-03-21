
################################################################################
# LIBRARIES #
################################################################################

library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(missMethyl)
library(writexl)
library(ggplot2)

################################################################################
# LOAD & ADAPT WHAT IS NEEDED TO ADAPT THE ANNOTATION #
################################################################################

## ANNOTATION
anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38::IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno$Name <- sub("_.*", "", anno$Name)

anno_red <- as.data.frame(unique(anno[, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "chr", "pos")]))
rownames(anno_red) <- c(1:nrow(anno_red))

# get rid of the non-CpG sites
anno_red<-anno_red[grepl("^cg",anno_red$Name),]


# REFGENE INFO FROM REFGENE (UCSC GENOME BROWSER, HG38)
genes <- read.table("OneDrive - University of Bergen/PhD/14_Downloads/Human_genes_RefSeq")
colnames(genes) <- c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
                     "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")

genes_red <- genes[, c("name2", "chrom", "strand", "txStart", "txEnd", "exonStarts", "exonEnds", "exonCount")]


################################################################################
# ADAPT THE ANNOTATION
################################################################################

for (i in 1:nrow(anno_red)) {
  
  if (i %% 10000 == 0) {
    print(i)
  }
  
  cg_anno <- as.character(anno_red$UCSC_RefGene_Name[i])
  
  if (cg_anno == "") {
    cg_name <- as.character(anno_red$Name[i])
    cg_pos <- as.numeric(anno_red$pos[i])
    cg_chr <- as.character(anno_red$chr[i])
    
    genes_filt <- genes_red[which(genes_red$chrom == cg_chr & genes_red$txStart <= cg_pos & genes_red$txEnd >= cg_pos), ]
    
    if (nrow(genes_filt) > 0) {
      gene_string <- paste0(genes_filt$name2, collapse = ";")
      anno_red[i, "UCSC_RefGene_Name"] <- gene_string
      
      gene_location_vector <- c()
      
      for (k in 1:nrow(genes_filt)) {
        exon_starts <- as.numeric(unlist(strsplit(genes_filt$exonStarts[k], ",")))
        exon_stops  <- as.numeric(unlist(strsplit(genes_filt$exonEnds[k], ",")))
        
        inside_exon <- which(cg_pos >= exon_starts & cg_pos <= exon_stops)
        
        if (cg_pos < exon_starts[1]) { 
          # 5'UTR
          gene_location_vector <- c(gene_location_vector, "5UTR")
        } else if (cg_pos > exon_stops[length(exon_stops)]) { 
          # 3'UTR
          gene_location_vector <- c(gene_location_vector, "3UTR")
        } else if (length(inside_exon) > 0) { 
          # Inside exon
          gene_location_vector <- c(gene_location_vector, paste0("exon_", inside_exon))
        } else { 
          # Intron
          intron_nr <- which(exon_starts > cg_pos)[1]
          gene_location_vector <- c(gene_location_vector, paste0("intron_", intron_nr))
        }
      }
      
      anno_red[i, "UCSC_RefGene_Group"] <- paste0(gene_location_vector, collapse = ";")
      
    } else {
      anno_red[i, "UCSC_RefGene_Name"] <- ""
      anno_red[i, "UCSC_RefGene_Group"] <- ""
    }
  }
}

anno_red$UCSC_RefGene_Group <- ifelse(anno_red$UCSC_RefGene_Name == "", "", anno_red$UCSC_RefGene_Group)

table(anno_red$UCSC_RefGene_Name == "")

#FALSE   TRUE 
#734946 185080 

save(anno_red, file = "OneDrive - University of Bergen/PhD/24_Case_Control_project/ANNOTATION_ADAPTED.RData")


################################################################################
# LOAD AND ADAPT SUMMARY STATS FOR GO #
################################################################################

load("OneDrive - University of Bergen/PhD/24_Case_Control_project/output/RData/DMPs_CaseCtrl_15CTRLPCs_2AncPCs_limma_anno_BACON.RData")
rownames(anno_red) <- anno_red$Name
all_cpgs <- DMPs_anno$cgID
ordered_DMPs_anno <- DMPs_anno[order(DMPs_anno$P.Value, decreasing = FALSE), ]
sign_0001 <- ordered_DMPs_anno$cgID[which(ordered_DMPs_anno$P.Value < 0.001)]

################################################################################
# RUN GO AND ADAPT RESULTS #
################################################################################
GO_result <- gometh(sig.cpg = sign_0001, all.cpg = all_cpgs, collection = "GO", array.type = "EPIC", genomic.features = "ALL", anno = anno_red)
GO_result$logP <- -log10(GO_result$P.DE)
GO_result <- GO_result[order(GO_result$P.DE), ]


################################################################################
# VISUALIZE & EXPORT RESULTS #
################################################################################

# import ggplot themes
load("OneDrive - University of Bergen/PhD/14_Downloads/themes/theme_transparent.RData")
load("OneDrive - University of Bergen/PhD/14_Downloads/themes/theme.RData")

# top 10 results
result_red <- GO_result[1:10, ]
result_red$significant <- ifelse(result_red$FDR <= 0.05, "Y", "N")

GO_plot <- ggplot(result_red, aes(x = logP, y = reorder(TERM, logP), fill = significant)) +
  geom_bar(stat = "identity") +
  labs(x = "-logP", y = NULL, title = NULL) +
  th + th_transparent +
  scale_fill_manual(values = c("Y" =  "#192f2d","N" = "#B0D0CE")) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0))
GO_plot

# export
ggsave("OneDrive - University of Bergen/PhD/24_Case_Control_project/GO_plot.pdf", plot = GO_plot, device = "pdf", width = 27, height = 10, units = "cm", dpi = 400)
write_xlsx(result_red, path = "OneDrive - University of Bergen/PhD/24_Case_Control_project/GO_table.xlsx")
