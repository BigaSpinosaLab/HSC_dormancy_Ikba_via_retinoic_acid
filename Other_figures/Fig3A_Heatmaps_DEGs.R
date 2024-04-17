################################################################################
##  Title       : RNAseq IkBa HSC dormancy - Thambyrajah R. et al. 2024
##                "IκBα controls dormancy in Hematopoietic Stem Cells via 
##                retinoic acid during embryonic development"
##  Description : This script is for plotting DEGs Heatmap (Fig3A)
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

require(openxlsx)
require(ComplexHeatmap)
require(dplyr)

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
                      
# Genes to highlight in the heatmap
genes <- c("Adgrg1","Gata2","Neurl3", "Sox18","Smad6","Rara", "Rarg", "Rarb", "Tfeb", "Cdk6",
           "Runx1","Gfi1","Itga2b") 
genes.ensembl <- final.res$gene_id[which(final.res$gene_name %in% genes)]

################################################################################
## Heatmap with DEGs
################################################################################

# FIRST, obtain Expression matrix with the genes of interest
mat <- norm.mat[which(rownames(norm.mat) %in% degs.sig$gene_id),] # These are the normalized counts for all samples

# Scale the rows (genes)
mat.z <- t(scale(t(mat))) # Equivalent to  t(apply(mat,1,scale))

# Value range in the matrix
quantile(mat.z, c(0.1, 0.95))
# 10%        95% 
# -1.142019  1.732627 # If abs(logFC) >1 is applied for selecting genes 

col_fun = circlize::colorRamp2(c(-1.5, 
                                 0, 
                                 2), c("#FF00FF", "black", "#FFFF00"))  # This is yellow/purple colors

cl = kmeans(mat.z, centers = 2)$cluster #Up and Down regulated in KO

htmp <- Heatmap(mat.z, name = "Scaled Exprs",  
                cluster_columns = TRUE,
                clustering_distance_columns = "euclidean", # Option by default
                show_column_dend = TRUE,
                
                #cluster_column_slices = TRUE,
                column_title_gp = gpar(fontsize = 10, fontface="bold"),
                column_gap = unit(0.5, "mm"),
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                show_row_names = FALSE,
                row_split =cl,
                row_title = c("Up-regulated", "Down-regulated"),
    
                col = col_fun,
                row_names_gp = gpar(fontsize = 9, fontface="bold"),
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

index =  which(rownames(mat.z) %in% genes.ensembl)
labels = rownames(mat.z)[index]

gnames <- final.res[which(final.res$gene_id %in% labels),c(1,20)]
gnames <- gnames[match(labels, gnames$gene_id),]  # Sort them in the same order as in labels

htmp <- htmp + rowAnnotation(sig_gene = anno_mark(at = index, labels = gnames$gene_name,
                                             side = "right", labels_gp = gpar(fontsize = 6, fontface="bold"),
                                             link_width = unit(3, "mm"),
                                             padding = unit(0, "mm"))) 

pdf("Fig3A_Heatmap_DEGs.pdf", width=3, height=5.5)
draw(htmp, heatmap_legend_side="bottom")
dev.off()
