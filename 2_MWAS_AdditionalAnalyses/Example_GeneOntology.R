
################################################################################
# LIBRARIES #
################################################################################

library(missMethyl)
library(writexl)
library(ggplot2)

################################################################################
# LOAD & ADAPT ANNOTATION
################################################################################

load("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/annotation/UPDATED_ANNOTATION.RData")
anno_upd <- unique(anno_upd[, c("Name", "Gene", "RegulRegion", "CHR", "MAPINFO")])
colnames(anno_upd) <- c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "chr", "pos")
anno_upd <- anno_upd[which(anno_upd$chr != "chr0"), ]
rownames(anno_upd) <- anno_upd$Name

################################################################################
# LOAD AND ADAPT SUMMARY STATS FOR GO #
################################################################################

load("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/DMPs/CaseCtrl_Meta/DMPs_meta_CaseCtrl_anno.RData")

all_cpgs <- DMPs_anno$cgID
ordered_DMPs_anno <- DMPs_anno[order(DMPs_anno$P.Value, decreasing = FALSE), ]
sign05percent <- round(nrow(DMPs_anno)*0.005, 0)

signDMPs <- ordered_DMPs_anno$cgID[1:sign05percent]

################################################################################
# RUN GO AND ADAPT RESULTS #
################################################################################

gometh <- function (sig.cpg, all.cpg = NULL, collection = c("GO", "KEGG"), 
          array.type = c("450K", "EPIC", "EPIC_V2"), plot.bias = FALSE, 
          prior.prob = TRUE, anno = NULL, equiv.cpg = TRUE, fract.counts = TRUE, 
          genomic.features = c("ALL", "TSS200", "TSS1500", "Body", 
                               "1stExon", "3'UTR", "5'UTR", "ExonBnd"), sig.genes = FALSE) 
{
  array.type <- match.arg(toupper(array.type), c("450K", "EPIC", 
                                                 "EPIC_V2"))
  collection <- match.arg(toupper(collection), c("GO", "KEGG"))
  genomic.features <- match.arg(genomic.features, c("ALL", 
                                                    "TSS200", "TSS1500", "Body", "1stExon", "3'UTR", "5'UTR", 
                                                    "ExonBnd"), several.ok = TRUE)
  if (length(genomic.features) > 1 & any(grepl("ALL", genomic.features))) {
    message("All input CpGs are used for testing.")
    genomic.features <- "ALL"
  }
  if (array.type == "450K" & any(grepl("ExonBnd", genomic.features))) {
    stop("'ExonBnd' is not an annotated feature on 450K arrays,\n\n           please remove it from your genomic.feature parameter\n\n           specification.")
  }
  if (array.type == "EPIC_V2" & any(grepl("ExonBnd", genomic.features))) {
    stop("'ExonBnd' is not an annotated feature on EPIC_V2 arrays,\n\n           please remove it from your genomic.feature parameter\n\n           specification.")
  }
  if (collection == "GO") {
    
      if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
        stop("org.Hs.eg.db package required but not installed.")
      egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
      GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
      d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
      GeneID.PathID <- GeneID.PathID[d, ]
      GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
                                                          keys=unique(GeneID.PathID$go_id), 
                                                          columns=c("GOID","ONTOLOGY","TERM"), 
                                                          keytype="GOID"))
      go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
      
      go <- list(idList=go, idTable=GOID.TERM)
      
      go$idTable <- go$idTable[which(go$idTable$ONTOLOGY == "BP"), ]
      
      go$idList <- go$idList[names(go$idList) %in% go$idTable$GOID]
      
    result <- gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = go$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(go$idTable, result, by.x = "GOID", by.y = "row.names")
    rownames(result) <- result$GOID
  }
  else if (collection == "KEGG") {
    kegg <- .getKEGG()
    result <- gsameth(sig.cpg = sig.cpg, all.cpg = all.cpg, 
                      collection = kegg$idList, array.type = array.type, 
                      plot.bias = plot.bias, prior.prob = prior.prob, anno = anno, 
                      equiv.cpg = equiv.cpg, fract.counts = fract.counts, 
                      genomic.features = genomic.features, sig.genes = sig.genes)
    result <- merge(kegg$idTable, result, by.x = "PathwayID", 
                    by.y = "row.names")
    rownames(result) <- result$PathwayID
  }
  result[, -1]
}



GO_result <- gometh(sig.cpg = signDMPs, all.cpg = all_cpgs, collection = "GO", array.type = "EPIC", genomic.features = "ALL", anno = anno_upd)
GO_result$logP <- -log10(GO_result$P.DE)
GO_result <- GO_result[order(GO_result$P.DE), ]


################################################################################
# VISUALIZE & EXPORT RESULTS #
################################################################################

# import ggplot themes
load("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/14_Downloads/themes/theme_transparent.RData")
load("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/14_Downloads/themes/theme.RData")

# top 5 results
result_red <- GO_result[1:5, ]
result_red$significant <- ifelse(result_red$FDR <= 0.05, "Y", "N")

GO_plot <- ggplot(result_red, aes(x = logP, y = reorder(TERM, logP), fill = significant)) +
  geom_bar(stat = "identity") +
  labs(x = "-logP", y = NULL, title = NULL) +
  th + th_transparent +
  ggtitle("GO Enrichment for the Top 0.5% Differentially Methylated CpGs") +
  scale_fill_manual(values = c("Y" =  "#192f2d","N" = "#B0D0CE")) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0))
GO_plot

# export
ggsave("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/DMPs/CaseCtrl_Meta/GO_plot.pdf", plot = GO_plot, device = "pdf", width = 25, height = 5, units = "cm", dpi = 400)
write_xlsx(result_red, path = "/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/DMPs/CaseCtrl_Meta/GO_table.xlsx")
