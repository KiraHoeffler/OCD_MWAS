
library(biomaRt)

#### ADAPT ANNOTATION ###
setwd("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/annotation/")

old_annotation <- read.csv("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/14_Downloads/MethylationEPIC v2.0 Files/EPIC-8v2-0_A1_manifest.csv", skip = 7)

ENCODE <- read.table("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/24_Case_Control_project/annotation/GRCh38-Closest-Genes-All.tsv", header = FALSE)
ENCODE <- ENCODE[, c(1:8, 10, 12, 13)]
colnames(ENCODE) <- c("chr", "start", "end", "cCRE_ID", "Experiment", 
                     "RegulElement", "TSS_chr", "TSS_pos", "transcriptID", "strand", "geneID")
ENCODE$genes_short <- sub("\\..*", "", ENCODE$geneID)




###############################################################################
                # TRANSLATE ENSEMBL IDs INTO GENE NAMES
################################################################################

write.table(ENCODE$genes_short, file = "IDs_before_translation.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)




#ADD GENE NAMES
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get the mappings
result <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "ensembl_gene_id",
                values = ENCODE$genes_short,
                mart = ensembl)


ensembl_grch37 <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  host = "https://grch37.ensembl.org"
)

# get gene names for old Ensembl IDs
your_ids <- c(ENCODE_incomplete$genes_short)

result_grch37 <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = your_ids,
  mart = ensembl_grch37
)


result_final <- rbind(result, result_grch37)

head(result_final)
save(result_final, file = "translation.RData")

ENCODE_merged <- merge.data.frame(x = ENCODE, y = result_final, by.x = "genes_short", by.y = "ensembl_gene_id", all.x = TRUE)
save(ENCODE_merged, file = "ENCODE_merged.RData")

result_final_empty <- result_final[which(result_final$external_gene_name == ""),]
ENCODE_empty <- ENCODE_merged[which(ENCODE_merged$external_gene_name == ""),]

writexl::write_xlsx(result_final_empty, path = "empty_translation.xlsx")

ENCODE_merged <- ENCODE_merged[which(ENCODE_merged$external_gene_name != ""), ]


ENCODE_merged$RegulElement2 <- NA
ENCODE_merged$RegulElement2[which(ENCODE_merged$RegulElement == "PLS")] <- "promoter"
ENCODE_merged$RegulElement2[which(ENCODE_merged$RegulElement == "pELS")] <- "prox_enhancer"
ENCODE_merged$RegulElement2[which(ENCODE_merged$RegulElement == "dELS")] <- "distal_enhancer"
ENCODE_merged$RegulElement2[which(ENCODE_merged$RegulElement == "CA-H3K4me3")] <- "promoter"

ENCODE_merged$RegulElement2[which(ENCODE_merged$RegulElement == "CA-TF")] <- "other_regul_region"
ENCODE_merged$RegulElement2[which(ENCODE_merged$RegulElement == "CA-CTCF")] <- "other_regul_region"
ENCODE_merged$RegulElement2[which(ENCODE_merged$RegulElement == "TF")] <- "other_regul_region"
ENCODE_merged$RegulElement2[which(ENCODE_merged$RegulElement == "CA")] <- "other_regul_region"


old_anno_red <- unique(old_annotation[, c("Name", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island")])
old_anno_red$Gene <- NA
old_anno_red$RegulRegion <- NA


for (i in 1:nrow(old_anno_red)){
  if (i %% 10000 == 0) {
    print(i)
  }
  
  chr <- old_anno_red$CHR[i]
  pos <- old_anno_red$MAPINFO[i]
  
  ENCODE_info <- unique(ENCODE_merged[which(ENCODE_merged$chr == chr & ENCODE_merged$start <= pos & ENCODE_merged$end >= pos), c("external_gene_name", "RegulElement2")])

  old_anno_red$Gene[i] <- paste0(ENCODE_info$external_gene_name, collapse = ";")
  old_anno_red$RegulRegion[i] <- paste0(ENCODE_info$RegulElement2, collapse = ";")
}



# REFGENE INFO FROM REFGENE (UCSC GENOME BROWSER, HG38)
genes <- read.table("/Users/kirahoeffler/OneDrive - University of Bergen/PhD/14_Downloads/Human_genes_RefSeq")
colnames(genes) <- c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
                     "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")

genes_red <- genes[, c("name2", "chrom", "strand", "txStart", "txEnd", "exonStarts", "exonEnds", "exonCount")]


################################################################################
# ADAPT THE ANNOTATION
################################################################################

for (i in 1:nrow(old_anno_red)) {
  
  if (i %% 10000 == 0) {
    print(i)
  }
  
  cg_anno <- as.character(old_anno_red$Gene[i])
  
  if (cg_anno == "") {
    cg_name <- as.character(old_anno_red$Name[i])
    cg_pos <- as.numeric(old_anno_red$MAPINFO[i])
    cg_chr <- as.character(old_anno_red$CHR[i])
    
    genes <- c()
    location <- c()
    
    #200
    genes_filt200plus <- genes_red[which(genes_red$chrom == cg_chr &
                                           genes_red$strand == "+" &
                                       cg_pos >= (genes_red$txStart - 200) & 
                                       cg_pos <= genes_red$txStart), ]
    
    genes_filt200minus <- genes_red[which(genes_red$chrom == cg_chr &
                                           genes_red$strand == "-" &
                                           cg_pos >= (genes_red$txEnd) & 
                                           cg_pos <= genes_red$txEnd + 200), ]
    genes_filt200 <- rbind(genes_filt200plus, genes_filt200minus)
    
    unique_genes <- unique(genes_filt200$name2)
    
    genes <- c(genes, unique_genes)
    location <- c(location, rep("TSS200", length(unique_genes)))
    
    #1500
    genes_filt1500plus <- genes_red[which(genes_red$chrom == cg_chr &
                                           genes_red$strand == "+" &
                                           cg_pos >= (genes_red$txStart - 1500) & 
                                           cg_pos <= genes_red$txStart), ]
    
    genes_filt1500minus <- genes_red[which(genes_red$chrom == cg_chr &
                                            genes_red$strand == "-" &
                                            cg_pos >= (genes_red$txEnd) & 
                                            cg_pos <= genes_red$txEnd + 1500), ]
    genes_filt1500 <- rbind(genes_filt1500plus, genes_filt1500minus)
    
    unique_genes2 <- unique(genes_filt1500$name2)
    
    genes <- c(genes, unique_genes2)
    location <- c(location, rep("TSS1500", length(unique_genes2)))
    
    
    # exon1
    genes_filt <- genes_red[which(genes_red$chrom == cg_chr & genes_red$txStart <= cg_pos & genes_red$txEnd >= cg_pos), ]
    genes_filt_plus <- genes_filt[which(genes_filt$strand == "+"), ]
    
    if (nrow(genes_filt_plus >= 1)){
      for (m in 1:nrow(genes_filt_plus)){
        exon_starts_plus <- as.numeric(unlist(strsplit(genes_filt_plus$exonStarts[m], ",")))
        exon_stops_plus  <- as.numeric(unlist(strsplit(genes_filt_plus$exonEnds[m], ",")))
        
        if (cg_pos <= exon_stops_plus[1] & cg_pos >= exon_starts_plus[1]){
          exon1_gene_plus <- genes_filt_plus$name2[m]
          genes <- c(genes, exon1_gene_plus)
          location <- c(location, "exon1")
        }
      }
    }
    
    genes_filt_minus <- genes_filt[which(genes_filt$strand == "-"), ]
    
    if (nrow(genes_filt_minus >= 1)){
      for (m in 1:nrow(genes_filt_minus)){
        exon_starts_minus <- as.numeric(unlist(strsplit(genes_filt_minus$exonStarts[m], ",")))
        exon_stops_minus  <- as.numeric(unlist(strsplit(genes_filt_minus$exonEnds[m], ",")))
        
        if (cg_pos <= exon_stops_minus[length(exon_stops_minus)] & cg_pos >= exon_starts_minus[length(exon_starts_minus)]){
          exon1_gene_minus <- genes_filt_minus$name2[m]
          genes <- c(genes, exon1_gene_minus)
          location <- c(location, "exon1")
        }
      }
    }
    
    gene_df <- unique(data.frame(genes = genes, location = location))
    old_anno_red$RegulRegion[i] <- paste0(gene_df$location, collapse = ";")
    old_anno_red$Gene[i] <- paste0(gene_df$genes, collapse = ";")
    
    
    if (old_anno_red$Gene[i] == ""){
      #gene body
      unique_genes <- unique(genes_filt$name2)
      
      old_anno_red$RegulRegion[i] <- paste0(rep("GeneBody", length(unique_genes)), collapse = ";")
      old_anno_red$Gene[i] <- paste0(unique_genes, collapse = ";")
    }
  }
}
   
old_anno_red$RegulRegion <- ifelse(old_anno_red$Gene == "", "", old_anno_red$RegulRegion)

table(old_anno_red$Gene == "")

anno_upd <- old_anno_red
save(anno_upd, file = "UPDATED_ANNOTATION.RData")

> table(old_anno_red$UCSC_RefGene_Name == "")

# FALSE   TRUE 
# 320874 609802 

