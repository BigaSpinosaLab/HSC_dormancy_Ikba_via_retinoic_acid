################################################################################
##  Title       : RNAseq IkBa HSC dormancy - Thambyrajah R. et al. 2024
##                "A critical requirement for IkBa in controlling dormancy in HSC
##                 via retinoic acid during embryonic development"
##  Description : This script is for conducting GSEA Pre ranked analysis against
##                general databases
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# NOTES:
# For conducting this analysis, it is required to have the results from Differential
# Expression Analysis (either execute 5_RNAseq_DEA.R script or use 'DEA_results
# _KO_vs_WT_HET' sheet from Suppl Table T3)

# REMARK: Obtained results are included in Suppl Table T3

################################################################################
## 1. Load packages
################################################################################

require(fgsea)
require(org.Mm.eg.db)
require(openxlsx)
require(msigdbr)
require(dplyr)
require(stringr)

################################################################################
## 2. Import data: complete list of genes with information from DEA
################################################################################

# Indicate path where DEA results are stored
# In this case, results enclosed in Suppl Table T3 are imported
path <- "../downloads/Suppl_Table_T4.xlsx"

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
## 3. Gene sets to be used: MSigDB Hallmark and ChEA TF
################################################################################

# MSigdB HallMark Gene sets
# Hallmark collection (others can be imported also from MSigdb)
# Check them: https://www.gsea-msigdb.org/gsea/msigdb/
msigdb_hallmark <- msigdbr(species = "Mus musculus", category = "H") 
# Put them in the required format
msigdbr_hallmark.gs = split(x = msigdb_hallmark$entrez_gene, f = msigdb_hallmark$gs_name)


# ChEA 2022 database
# Downloaded from: https://maayanlab.cloud/Enrichr/#libraries
chea22 <- GSA::GSA.read.gmt(filename = "../../AGM_Notch1/RNAseq_AGM_Notch1/R/ChEA_2022.txt")

# For GSEA, it must be a list of vectors, each vector should refer to a particular pathway
chea22.gsea <- chea22$genesets
names(chea22.gsea) <- chea22$geneset.names

# Select only Mouse datasets (there are human and mouse experiments)
chea22.gsea <- chea22.gsea[grep("Mouse", names(chea22.gsea))]

# In all gene sets there is an empty element in last position. Let's remove it
chea22.gsea <- lapply(chea22.gsea, function(gs) gs[-length(gs)])

# And finally, transform gene symbols into Entrez IDs

chea22.gsea <- lapply(chea22.gsea, function (gs) {
  transform <- clusterProfiler::bitr(geneID = str_to_title(gs), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop=TRUE)
  return(transform$ENTREZID)
})

# Store gene sets from different databases into a list
databases.gs <- list("Hallmark" = msigdbr_hallmark.gs,
                     "ChEA2022" = chea22.gsea)

################################################################################
## 4. Compute Pre-Ranked GSEA for KEGG gene sets
################################################################################

ranking <- ranked.genes.red 
  
# Compute fgsea for each database

fgseaRes.all <- lapply(databases.gs, function(database){
    set.seed(123)
    fgseaRes <- fgseaMultilevel(pathways = database, 
                                stats = ranking,
                                scoreType = "std",
                                minSize=10,
                                nPermSimple = 10000,
                                maxSize=500)
    
    # LeadingEdge column: Transform entrez ids into symbols
    fgseaRes[, leadingEdge2 := mapIdsList(
      x=org.Mm.eg.db, 
      keys=leadingEdge,
      keytype="ENTREZID", 
      column="SYMBOL")]
    
    # Return results
    return(fgseaRes)
})
  

## Store results in an excel and RData file
# Fix leading edge columns
fgseaRes.all <- lapply(fgseaRes.all, function(database){
    database <- database %>% 
      mutate(leadingEdge = sapply(leadingEdge, toString),
             leadingEdge2 = sapply(leadingEdge2, toString))
})
  
fgseaRes.all <- lapply(fgseaRes.all, function(res) res[order(res$padj),])
  
# Store resuls in an excel
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(x = fgseaRes.all, file= paste(path,"/GSEA_", suffix,".xlsx",sep=""),
             headerStyle=hs)
  
