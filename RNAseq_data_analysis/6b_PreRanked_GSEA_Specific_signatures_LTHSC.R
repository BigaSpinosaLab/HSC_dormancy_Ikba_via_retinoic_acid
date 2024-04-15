################################################################################
##  Title       : RNAseq IkBa HSC dormancy - Thambyrajah R. et al. 2024
##                "A critical requirement for IkBa in controlling dormancy in HSC
##                 via retinoic acid during embryonic development"
##  Description : This script is for GSEA Pre ranked analysis for particular gene
##                signatures: one from Manesia et al 2017 and other from 
##                Thambyrajah et al 2023 (biorxiv) or 2024 Nat Comm
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# NOTES:
# For conducting this analysis, it is required to have the results from Differential
# Expression Analysis (either execute 5_RNAseq_DEA.R script or use 'DEA_results
# _KO_vs_WT_HET' sheet from Suppl Table T3)

# REMARK: Obtained results are included in Suppl Table T3

# Gene set to be taken from Manesia et al 2017  (Fetal liver HSC)
# http://dx.doi.org/10.1089/scd.2016.0294

# Gene set to be taken from Thambyrajah et al 2023 (AGM HSC-HE /HSC)
# https://www.biorxiv.org/content/10.1101/2023.04.19.537430v1

################################################################################
## 1. Load packages
################################################################################

require(fgsea)
require(KEGGREST)  # To be used for retrieving the KEGG and GO gene sets
require(org.Mm.eg.db)
require(clusterProfiler)
require(openxlsx)
require(dplyr)
require(grid)  # For adding some annotations in the plot
require(data.table) # For the adapted plot

require(ggvenn)

################################################################################
## 2. Import data: complete list of genes 
################################################################################

# Indicate path where DEA results are stored
# In this case, results enclosed in Suppl Table T3 are imported
path <- "../downloads/Suppl_Table_T3.xlsx"

# DEGs obtained from comparing KO vs (HET + WT) - Covariate adjustment (Gender)
final.res <- read.xlsx(xlsxFile = path, sheet = "DEA_results_KO_vs_WT_HET")
suffix="KO_vs_HETWT" #For naming purposes

# Rank genes based on log2 FC shrunken and using the ensembl gene ids
ranked.genes <- final.res$Shrunken_lFC
names(ranked.genes) <- final.res$Entrez
ranked.genes <- sort(ranked.genes, decreasing=TRUE) 

# Genes with no Entrez annotation  out of the initial 20,610 genes
length(which(is.na(names(ranked.genes))))
# [1] 4,190

# Remove those ranked genes without any name
ranked.genes.red <-  ranked.genes[-which(is.na(names(ranked.genes)))]

# Remove duplicated gene names. Check that these are exceptionally (too few values)
length(which(duplicated(names(ranked.genes.red))))  # 5
ranked.genes.red <- ranked.genes.red[-which(duplicated(names(ranked.genes.red)))]

################################################################################
## 4. Import data: Gene sets to be used
################################################################################

#===============
# GeneSet A: Fetal Liver LT-HSC signature from Manesia et al. 2017
#===============
manesia_etal.degs <- read.xlsx(xlsxFile = "Manesia_etal_2017_DEGs_FL_BM.xlsx")  # All genes

manesia_etal.degs <- manesia_etal.degs %>%
  filter(P.adj <0.05, abs(log2FoldChange) >1) # Criteria used by Manesia et al. 

# Take only the top300 (sorted by lFC) DEGs with logFC <0 (UPREGULATED in FL)
# interpreted as the signature 

manesia.sig <- manesia_etal.degs %>% 
  filter(log2FoldChange<0) %>%
  arrange(log2FoldChange)

manesia.sig <- manesia.sig[1:300,"ENTREZID"]  

#===============
# GeneSet B: AGM HSC-HE and HSC signature from Thambyrajah et al. 2023 (cluster representatives from scRNAseq data)
# or Thambyrajah et al. 2024 
#===============
thambyrajah_etal.repHSCHE <- read.xlsx(xlsxFile = "Thambyrajah_etal_2024_DEGs_all_clusters.xlsx", 
                                   sheet="Gene_marker_HSC-HE")
thambyrajah_etal.repHSC <- read.xlsx(xlsxFile = "Thambyrajah_etal_2024_DEGs_all_clusters.xlsx", 
                                       sheet="Gene_marker_HSC")

# First, reduce the list of genes to those with: cells percentage in the specific cluster > 50%;
# lFC >0 and adj pval (<1e-5)

thambyrajah_etal.repHSCHE.sig <- thambyrajah_etal.repHSCHE %>%
  filter(HE.pct > 0.5, adj.pval < 0.00001, HE_l>0)  %>%
  arrange(desc(HE_l))  

thambyrajah_etal.repHSCHE.sig <- clusterProfiler::bitr(geneID = thambyrajah_etal.repHSCHE.sig$HE_n, fromType = "SYMBOL", 
                                                       toType = "ENTREZID", 
                                                       OrgDb = org.Mm.eg.db, 
                                                       drop = TRUE)

thambyrajah_etal.repHSCHE.sig <- thambyrajah_etal.repHSCHE.sig$ENTREZID[1:150]


thambyrajah_etal.repHSC.sig <- thambyrajah_etal.repHSC %>%
  filter(HSC.pct > 0.5, adj.pval < 0.00001, HSC_l>0)  %>%
  arrange(desc(HSC_l))  
thambyrajah_etal.repHSC.sig <- clusterProfiler::bitr(geneID = thambyrajah_etal.repHSC.sig$HSC_n, fromType = "SYMBOL", 
                                                       toType = "ENTREZID", 
                                                       OrgDb = org.Mm.eg.db, 
                                                       drop = TRUE)

thambyrajah_etal.repHSC.sig <- thambyrajah_etal.repHSC.sig$ENTREZID[1:150]

# Combine both and transfrom gene names into Entrez
thambyrajah_etal.rep <- unique(c(thambyrajah_etal.repHSCHE.sig, thambyrajah_etal.repHSC.sig))

#===============
# GeneSet C: NfKB targets (current Suppl Fig4) + Il6 and Rante (the two latter are not present in our matrix)
#===============
genes <- c("Tnfaip3","Irf1","Irf2", "Myb","Myc","Dctn3", "Sqstm1", "Ripk1", "Il1a", "Il1b","Il1r2","Tradd", "Il6","Ccl5")
genes.ensembl <- clusterProfiler::bitr(geneID = genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=TRUE)

#===============
# GeneSet D: genes from KEGG Pathway mmu04064 (NFkb signaling pathway)
#===============
names <- keggGet("mmu04064")[[1]]$GENE
# Entrez ids
nfkb.entrez <- names[seq(1,length(names),2)]

#===============
# Put all gene sets into one object
#===============

custom.gs <- list("Manesia_FL_LTHSC" = manesia.sig,
                  "Thambyrajah_AGM_HSC_HE" = thambyrajah_etal.rep,
                  "NfkB_custom_targets" = genes.ensembl$ENTREZID,
                  "KEGG_NfkB_signaling_pathway" = nfkb.entrez)


################################################################################
## Explore the overlapping between these two signatures (FL and AGM)
################################################################################

tovenn <- custom.gs
names(tovenn) <- c("FL LT-HSC", "AGM HSC-HE")

pdf(file="Venn_Custom_GS_FL_LTHSC_AGM_HSC_common_Entrez.pdf",width=1.5,height=1.5)
ggvenn(
  tovenn, 
  fill_color = c('brown2', 'green3'),
  set_name_size = 2,
  text_size = 1.8,
  stroke_color = "white",stroke_size = 0.1,
  fill_alpha = 0.3) 
dev.off()

################################################################################
## Compute Pre-Ranked GSEA for custom gene sets
################################################################################

set.seed(123)
fgseaRes <- fgseaMultilevel(pathways = custom.gs, 
                            stats = ranked.genes.red ,
                            scoreType = "std",
                            nPermSimple = 10000)
    
# LeadingEdge column: Transform entrez ids into symbols
fgseaRes[, leadingEdge2 := mapIdsList(
      x=org.Mm.eg.db, 
      keys=leadingEdge,
      keytype="ENTREZID", 
      column="SYMBOL")]
    
## Store results in an excel and RData
# Fix leading edge columns
fgseaRes <- fgseaRes %>% 
      mutate(leadingEdge = sapply(leadingEdge, toString),
             leadingEdge2 = sapply(leadingEdge2, toString))

# Store resuls in an excel
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(x = fgseaRes, 
           file= "GSEA_Random_Walk_Custom_Signatures.xlsx",
             headerStyle=hs)

################################################################################
## Enrichment plot with results
################################################################################

pdf(file = "GSEA_plot_FL_LTHSC_Magnesia.pdf", width=4, height=4)
plotEnrichment_adapted(pathway = custom.gs$Manesia_FL_LTHSC, 
                       stats = ranked.genes.red,
                       col_walking = "brown2",
                       fgsea_res = fgseaRes[1,],
                       ypos=0.38) + 
  labs(title="Fetal Liver LT-HSC (Manesia et al.)") 
dev.off()


pdf(file = "GSEA_plot_AGM_Thambyrajah.pdf", width=4, height=4)
plotEnrichment_adapted(pathway = custom.gs$Thambyrajah_AGM_HSC_HE, 
                       stats = ranked.genes.red,
                       fgsea_res = fgseaRes[3,],
                       col_walking = "green3",
                       ypos = 0.2) + 
  labs(title="AGM HSC-HE (Thambyrajah et al.)") 
dev.off()


pdf(file = "GSEA_plot_NfkB_targets_selected.pdf", width=4, height=4)
plotEnrichment_adapted(pathway = custom.gs$NfkB_custom_targets, 
                       stats = ranked.genes.red,
                       fgsea_res = fgseaRes[2,],
                       col_walking = "blue",
                       ypos = 0.2) + 
  labs(title="NfKB targets - selected") 
dev.off()

pdf(file = "GSEA_plot_NfkB_signaling_pathway_KEGG.pdf", width=4, height=4)
plotEnrichment_adapted(pathway = custom.gs$KEGG_NfkB_signaling_pathway, 
                       stats = ranked.genes.red,
                       fgsea_res = fgseaRes[1,],
                       col_walking = "blue",
                       ypos = 0.2) + 
  labs(title="NfKB signaling (mmu04064)") 
dev.off()

################################################################################
## ANNEX. Enrichment plot adapted (aesthetics) from the original from fgsea package
################################################################################

require(data.table)  

plotEnrichmentData <- function(pathway, stats,
                                 gseaParam=1) {
    
    if (any(!is.finite(stats))){
      stop("Not all stats values are finite numbers")
    }
    
    rnk <- rank(-stats)
    ord <- order(rnk)
    
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    pathway <- unique(pathway)
    
    gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)
    
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.table(rank=c(0, xs, n + 1), ES=c(0, ys, 0))
    ticks <- data.table(rank=pathway, stat=statsAdj[pathway])
    stats <- data.table(rank=seq_along(stats), stat=statsAdj)
    
    res <- list(
      curve=toPlot,
      ticks=ticks,
      stats=stats,
      posES=max(tops),
      negES=min(bottoms),
      spreadES=max(tops)-min(bottoms),
      maxAbsStat=max(abs(statsAdj)))
  }
  

plotEnrichment_adapted <- function(pathway, 
                                   stats,
                                   fgsea_res,
                                   col_walking,
                                   ypos,
                             gseaParam=1,
                             ticksSize=0.2) {
    
    pd <- plotEnrichmentData(
      pathway = pathway,
      stats = stats,
      gseaParam = gseaParam)
    
    # Create the text to be added in the plot regarding NES and padj
    nes = round(fgsea_res$NES, digits=2)
    padj = formatC(fgsea_res$pval, digits=2, format="e")  # Let's take the pval since we have only tested two gene sets
    # It does not make sense to make FDR in this scenario
    
    grob <- grobTree(textGrob(paste("NES=", nes,"\npval=",padj,sep=""), x=0.6,  y=0.85, hjust=0,
                              gp=gpar(col="gray43", 
                                      fontsize=13, 
                                      fontface="bold")))

    # Create the text to include KO and HE-WT labels
    grob2 <- grobTree(textGrob("KO", x=0.1,  y=ypos, hjust=0,
                              gp=gpar(col="darkorchid4", 
                                      fontsize=15, 
                                      fontface="bold")))
    
    grob3 <- grobTree(textGrob("HE-WT", x=0.7,  y=ypos, hjust=0,
                               gp=gpar(col="khaki4", 
                                       fontsize=15, 
                                       fontface="bold")))
    
    
    with(pd,
         ggplot(data=curve) +
           geom_line(aes(x=rank, y=ES), color=col_walking) +  # We plot the NES instead of ES
           geom_segment(data=ticks,
                        mapping=aes(x=rank, y=-spreadES/16,
                                    xend=rank, yend=spreadES/16),
                        linewidth=ticksSize) +
           annotation_custom(grob) +
           annotation_custom(grob2) +
           annotation_custom(grob3) +
           geom_hline(yintercept=posES, colour="red", linetype="dashed") +
           geom_hline(yintercept=negES, colour="red", linetype="dashed") +
           geom_hline(yintercept=0, colour="black") +
           theme(
             #panel.border=element_rect(color="black"),
             panel.background = element_blank(),
             panel.grid.major =element_blank()
             #panel.grid.major=element_line(color="grey92")
           ) +
           labs(x="Rank", y="Enrichment Score"))
  }
  

