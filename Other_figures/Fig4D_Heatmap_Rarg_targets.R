################################################################################
##  Title       : RNAseq IkBa HSC dormancy - Thambyrajah R. et al. 2024
##                "IκBα controls dormancy in Hematopoietic Stem Cells via 
##                retinoic acid during embryonic development"
##  Description : This script is for plotting Rarg targets expression levels in 
##                a heatmap
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# REMARK: Rarg target genes listed in TRANSFAC database (id TF:M10353) are required
# These can be achieved from gprofiler2 (function gconvert)

# Supp Data T3 is required for importing data

require(openxlsx)
require(ComplexHeatmap)
require(dplyr)
require(gprofiler2)
require(org.Mm.eg.db)

################################################################################
## Import data and identify DEGs
################################################################################

# DEGs obtained from comparing KO vs (HET + WT) - Covariate adjustment (Gender)
# This can be obtained from Supp Data T3
final.res <- read.xlsx(xlsxFile = "Supp_Data_T3.xlsx", sheet="DEA_results_KO_vs_WT-HET",
                       startRow = 4) # To avoid header

# Identify DEGs per adj pval and shrunken logFC criteria
degs.sig <- final.res[which(final.res$padj<0.05),]
degs.sig <- degs.sig[which(abs(degs.sig$Shrunken_lFC)>1),]

# Import expression matrix (normalized and gender corrected)
norm.mat <- read.xlsx(xlsxFile = "Supp_Data_T3.xlsx", sheet="Exprs_mat_Norm_Gender_corrected")
rownames(norm.mat) <- norm.mat$gene_id
norm.mat <- norm.mat[,-1]

#colnames(norm.mat)
# HE1 HE2 HE3 KO1 KO3 WT1 WT2 WT3

# Metadata file
targets <- data.frame("sampleID" = colnames(norm.mat),
                      "condition" = c(rep("HE",3), rep("KO",2), rep("WT",3)),
                      "gender" = c("M", "F", "F", "F", "F", "M", "M", "F"))


################################################################################
## Import data: Target genes from TFs of interest (TRANSFAC)
################################################################################

# TFs of interest with their TRANSFAC ids
tfs <- list("Rarg" = "TF:M10353")

tfs.targets <- lapply(tfs, function(tf) gconvert(query = tf, organism = "mmusculus", target="ENSG"))
tfs.targets.ensembl <- lapply(tfs.targets, function(tf) tf$target)

# Alternatively, they can be imported thru xlsx file in repo

################################################################################
## Preamble: Transform ENSEMBL IDs in the expression matrix into Entrez
################################################################################

# Select which TF
experiment <- "Rarg"
genes.ensembl <- tfs.targets[[experiment]]

  # Remove duplicates by gene name
genes.ensembl <- genes.ensembl[!duplicated(genes.ensembl$name),]

# Subset expression matrix for the genes of interest
mat <- norm.mat[which(rownames(norm.mat) %in% genes.ensembl$target),] 

  # Rarg: 2868 out of 4326 included

# Use gene names instead of ENSMBL ids for row names -> to show names in the heatmap
genes.ensembl <- genes.ensembl[match(rownames(mat), genes.ensembl$target),]
rownames(mat) <- genes.ensembl$name

################################################################################
## Vertical panel indicating if those genes are DEGs or not
################################################################################

degs.to.plot <- data.frame("gene_id" = rownames(mat), 
                           "DEG" = "No")  

degs.to.plot$DEG[which(degs.to.plot$gene_id %in% degs.sig$gene_name)] <- "Yes" 

# Lateral annotation with DEGs
col_DEGs = c("Yes" = "darkolivegreen2", "No" = "dodgerblue4")
ht.degs <- Heatmap(degs.to.plot$DEG,
                      name ="DEG",
                      col = col_DEGs, 
                      width = unit(4, "mm"),
                      show_heatmap_legend = TRUE)

################################################################################
## Heatmap with genes of interest
################################################################################

# Scale the rows (genes)
mat.z <- t(scale(t(mat))) # Equivalent to  t(apply(mat,1,scale))

# Value range in the matrix
quantile(mat.z, c(0.1, 0.95))
# 10%        95% 
#-1.260844  1.563227  #Rara TF

col_fun = circlize::colorRamp2(c(-1.5, 0, 2), c("#FF00FF", "black", "#FFFF00"))  # This is yellow/purple colors

htmp <- Heatmap(mat.z, name = "Scaled Exprs", 
                column_title = paste0("TRANSFAC TF: ",experiment),
                cluster_columns = TRUE,
                clustering_distance_columns = "euclidean", # Option by default
                show_column_dend = TRUE,
                #cluster_column_slices = TRUE,
                column_title_gp = gpar(fontsize = 10, fontface="bold"),
                column_gap = unit(0.5, "mm"),
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                show_row_names = FALSE,
                row_names_side = "left",
                col = col_fun,
                row_names_gp = gpar(fontsize = 6, fontface="bold"),
                top_annotation = HeatmapAnnotation(Phenotype = factor(targets$condition),
                                                   Gender = factor(targets$gender),
                                                   col = list(Phenotype = c("HE" = "khaki4",
                                                                                  "WT" = "gray50",
                                                                                  "KO" = "darkorchid4"),
                                                              Gender = c("F" = "salmon2",
                                                                            "M" = "olivedrab3")),
                                                   simple_anno_size=unit(0.3,"cm"),
                                                   show_annotation_name = c(Phenotype=FALSE, Gender=FALSE)), 
                
                show_column_names = TRUE,
                use_raster = TRUE,
                heatmap_legend_param = list(legend_direction = "horizontal"),
                raster_quality = 4)

# To show the names of those genes that are DEGs
index =  which(rownames(mat.z) %in% degs.to.plot[degs.to.plot$DEG %in% "Yes","gene_id"])
labels = rownames(mat.z)[index]

ht.gnames <- rowAnnotation(sig_gene = anno_mark(at = index, labels = labels, #gnames$gene_name,
                                                  side = "right", 
                                                  labels_gp = gpar(fontsize = 3, fontface="bold"),  #It was 3 for the other TFs
                                                  link_width = unit(3, "mm"),
                                                  padding = unit(0, "mm"))) 

ht <- htmp + ht.degs + ht.gnames

pdf(paste0("Fig_4D_Heatmap_TRANSFAC_", experiment,".pdf"), width=3.5, height=11)  
draw(ht, heatmap_legend_side="bottom")
dev.off()

