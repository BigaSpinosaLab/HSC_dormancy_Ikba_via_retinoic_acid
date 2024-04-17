################################################################################
##        Title : RNAseq IkBa HSC dormancy - Thambyrajah R. et al. 2024
##                "IκBα controls dormancy in Hematopoietic Stem Cells via 
##                retinoic acid during embryonic development"
##  Description : This script is for conducting Differential Expression Analysis 
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# NOTES: 
# a) Required input data for this analysis is an expression matrix (RAW counts).
#    This matrix can be directly downloaded from GEO repository (GSE188523) or
#    created/reproduced using original FASTQ files + shared shell scripts. 
# b) Fitted statistical model with two covariates: gender and genotype considering
#    WT and HETs in the same group and as the reference. Thus, comparison is
#    KO vs WT-HET.
# c) By executing this script, DEA results in Suppl Table T3 are obtained
# d) Wrap-up function to gather DEA results is in the last section of this script
# (names as complete.dea())

# REMARK: Basic exploratory analysis is performed (PCA)
# REMARK: This script also includes the statistical analysis for comparing between
# WT and HET.

################################################################################
## 1. Load packages
################################################################################

require(DESeq2)
require(ggplot2)
require(openxlsx)
require(EnhancedVolcano)
require(vsn)
require(limma) # To correct for gender

################################################################################
## 2. Import relevant data and prepare it for analysis
################################################################################

# Raw count matrix from GEO repo > GSE188523_RNAseq_All_samples_Raw_counts.txt
# Gene quantification is directly obtained from STAR

counts <- read.delim(file="/Users/mmaqueda/Downloads/GSE188523_RNAseq_All_samples_Raw_counts.txt", 
                     header=TRUE, sep="\t")

# Create dataframe with samples metadata information (Gender and Genotype)
# colnames(counts)
# [1] "GeneId" "HE1"    "HE2"    "HE3"    "KO1"    "KO2"    "WT1"    "WT2"    "WT3"   

  # Define genes as row names
rownames(counts) <- counts$GeneId
counts <- counts[,-1]

  # Metadata file
metadata <- data.frame("Sample.id" = colnames(counts),
                      "Gender" = c("M", "F", "F", "F", "F", "M", "M", "F"),
                      "Condition" = c(rep("HET",3), 
                                      rep("KO",2), 
                                      rep("WT",3)))

# Add a join condition for HET and WT for further use
metadata$ReCondition <- ifelse(metadata$Condition == "KO", "KO", "HE_WT")


# Gene annotation (dim 55487x14): Created from corresponding GTF file from Ensembl
# Shown an entry example
gene.annot <- readRDS(file = "/Volumes/projectscomput/cancer/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.102.ONLY_GENES.RDS")

# chrom  source feature   start     end score strand frame                                                                                               attribute
# 1     ensembl    gene 3102016 3102125     .      +     . gene_id ENSMUSG00000064842; gene_version 1; gene_name Gm26206; gene_source ensembl; gene_biotype snRNA;

# gene_id gene_type gene_name    Entrez gene_biotype
# ENSMUSG00000064842      <NA>   Gm26206 115487594        snRNA

################################################################################
## 3. Create DESeq2 object with all samples to make exploratory analysis
## Remove very low expressed genes and apply VST norm exprs matrix
################################################################################

# Create object
data.dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                           colData = metadata, 
                                           design =  ~ Gender + Condition,
                                           tidy=FALSE)

dim(data.dds) 
#[1] 55487     8
colData(data.dds)
head(assay(data.dds,"counts"))

# Remove low-count genes before continue

  # Remove those with less than 10 counts accross all samples
data.dds.red <- data.dds[rowSums(counts(data.dds)) >=10]
nrow(data.dds.red)   
# [1] 20610 

  # Discard those annotations related to those genes
gene.annot.red <- gene.annot[rowSums(counts(data.dds)) >=10,]

  # Transformation: for visualization and clustering purposes

#colSums(assay(data.dds.red))
# HE1      HE2      HE3      KO1      KO2      WT1      WT2      WT3 
# 19390262 21406146 19711324 17891738 19940715 19605203 20635376 21716263 

vstData.dds <- vst(data.dds.red, blind=FALSE) 

################################################################################
## 4. Exploratory analysis: PCA with corrected gender
################################################################################

# Manual computation with prcomp
# ===

mat <- assay(vstData.dds) 
cols = c("WT" = "gray50", "HE" = "khaki4", "KO" = "darkorchid4")

# Let's plot PCA after removing the effect of Gender
mm <- model.matrix(~ condition, colData(vstData.dds))
mat2 <- limma::removeBatchEffect(x = mat, 
                                 batch=vstData.dds$gender,
                                 design=mm)
pca <-  prcomp(t(mat2), scale=TRUE, center=TRUE)
PCs <- as.data.frame(pca$x)
PCs$Sample <- rownames(PCs)
PCs$Condition <- as.character(colData(vstData.dds)$Condition)
PCs$Gender <- as.character(colData(vstData.dds)$Gender)
Variance <- round(summary(pca)$importance[2,]*100, digits=1)


pdf(file = "RNAseq_IkBa_HSC_PCA_Complete_Dataset_Gender_corrected.pdf", width=6, height=6)
ggplot(PCs, aes(PC1, PC2, fill=Condition,label=Sample, col= Condition, shape=Gender)) +
  geom_point(size=7,alpha=0.8,) +
  geom_text(aes(label=Sample),vjust= 1.8,hjust=0.5,size=3, color="black")+
  xlab(paste0("PC1: ",Variance[1],"% variance")) +
  ylab(paste0("PC2: ",Variance[2],"% variance")) +
  scale_color_manual(values = cols) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color="black"),
        legend.title = element_text(size = 12,face="italic"),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1)) +
  coord_fixed() +
  xlim(-150,200) +
  ylim(-100,100) +
  ggtitle("RNAseq HSC IkBa - Complete Dataset - Corrected gender") 
dev.off()

################################################################################
## 5. (Suppl) Are there significant differences between WT vs HETs? Conduct DEA
## to assess if WT and HETs can be combined
## SKIP THIS SECTION IF ONLY INTERESTED in KO vs WT-HETs comparison
################################################################################

# Create a new object only considering WT and HETs samples
data.dds.wt.hets <- DESeq2::DESeqDataSetFromMatrix(countData = counts[,which(colnames(counts) %in% c("HE1","HE2","HE3", "WT1","WT2","WT3"))], 
                                           colData = metadata[which(metadata$Sample.id %in% c("HE1","HE2","HE3", "WT1","WT2","WT3")),], 
                                           design =  ~ Gender + Condition,
                                           tidy=FALSE)

# Remove low-count genes before continue
  # Remove those with less than 10 counts accross all samples
data.dds.wt.hets.red <- data.dds.wt.hets[rowSums(counts(data.dds.wt.hets)) >=10]
nrow(data.dds.wt.hets.red)   
# [1] 19542

# Discard those annotations related to those genes
gene.annot.wt.hets.red <- gene.annot[rowSums(counts(data.dds.wt.hets)) >=10,]

# VST computation on data
vstData.dds.wt.hets <- vst(data.dds.wt.hets.red, blind=FALSE) 

# Relevel the condit
colData(data.dds.wt.hets.red)$Condition <- relevel(colData(data.dds.wt.hets.red)$Condition, ref = "WT")
colData(data.dds.wt.hets.red)$Gender <- relevel(colData(data.dds.wt.hets.red)$Gender, ref = "F")

# Apply DESeq2
dea.dds <- DESeq(data.dds.wt.hets.red)

# Get results
resultsNames(dea.dds)  #Identify the name of the comparison of interest
# [1] "Intercept"          "gender_M_vs_F"      "condition_HE_vs_WT"

comparisons <- c("condition_HE_vs_WT")  # Just one comparison

all.comparisons <- lapply(comparisons, function(comp){
  dea <- complete.dea(DESeqSet = dea.dds,
                      comparison = comp,
                      fc_threshold  = 1, 
                      gene_annotations = gene.annot.wt.hets.red,
                      suffix.fig.files = "_HE_vs_WT_subset",
                      lim.ma.abs = 10,
                      norm.data=vstData.dds.wt.hets,
                      res.path = "PAPER/Hets_vs_WT_subset_ANALYSIS/")
  return(dea)
})

names(all.comparisons) <- comparisons

################################################################################
## 6. FINAL DEA: aggregating HETs and WTs  - Complete dataset
################################################################################

# Create object again, now with other design
data.dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                           colData = metadata, 
                                           design =  ~ Gender + ReCondition,
                                           tidy=FALSE)

dim(data.dds) 
#[1] 55487     8

# Remove those with less than 10 counts accross all samples
data.dds.red <- data.dds[rowSums(counts(data.dds)) >=10]
nrow(data.dds.red)   
# [1] 20610 

# Discard those annotations related to those genes
gene.annot.red <- gene.annot[rowSums(counts(data.dds)) >=10,]

# Normalize data
vstData.dds <- vst(data.dds.red, blind=FALSE) 

# Relevel the condit
colData(data.dds.red)$ReCondition <- relevel(colData(data.dds.red)$ReCondition, ref = "HE_WT")
colData(data.dds.red)$Gender <- relevel(colData(data.dds.red)$Gender, ref = "F")

# Apply DESeq2
dea.dds <- DESeq(data.dds.red)
comparisons <- c("ReCondition_KO_vs_HE_WT")  

all.comparisons <- lapply(comparisons, function(comp){
  dea <- complete.dea(DESeqSet = dea.dds,
                      comparison = comp,
                      fc_threshold  = 1, 
                      gene_annotations = gene.annot.red,
                      suffix.fig.files = "_complete_dataset",
                      lim.ma.abs = 10,
                      norm.data=vstData.dds,
                      samples=c("HE","KO","WT"),
                      res.path = "PAPER/DEA_complete_dataset/")
  return(dea)
})

names(all.comparisons) <- comparisons

################################################################################
## ANNEX Function to compute DEA, lfcShrinkage, gather all information and add annotations
##    into the same data frame and generate a Volcano plot + excel file with results
################################################################################

complete.dea <- function(DESeqSet, #Obtain after computing DESeq function
                         comparison, 
                         gene_annotations, # File with the corresponding gene annotations 
                         suffix.fig.files,
                         fc_threshold, # Shrunken fold change threshold for Volcano Plot
                         lim.ma.abs= 10,  #Limit to consider in MA plot
                         res.path,  #Path to store MA plot and excel file with the results from DEA
                         samples = NULL, # If not provided, from comparison names (only valid for not aggregated comparison)
                         norm.data) # VST normalized data 
{
  # 1. Obtain results from DEA
  results_dea <-results(DESeqSet,
                        name =  comparison,
                        alpha = 0.05)
  
  # 2. Print summary of results on the console
  print(comparison)
  summary(results_dea)
  
  # 3. LogFC shrinkage
  res.shr <- lfcShrink(dds = DESeqSet, 
                       coef= comparison,
                       res = results_dea, 
                       type = "apeglm")
  
  # 4. Store an MA-plot to see the effect
  pdf(file = file.path(res.path,paste("MAplot_",comparison,suffix.fig.files,".pdf",sep="")), width=16, height=10)
  par(mfrow=c(1,2))    
  DESeq2::plotMA(results_dea,ylim=c(-lim.ma.abs,lim.ma.abs),main=paste("DEA ",comparison,suffix.fig.files, sep=""))
  DESeq2::plotMA(res.shr, ylim=c(-lim.ma.abs,lim.ma.abs), main=paste("DEA ", comparison,suffix.fig.files," - Shrinkage performed", sep=""))
  dev.off()
  
  # 5. Collect results in the same data frame: DEA resuls and logFC shrinkage
  DEA_complete <- cbind(as.data.frame(results_dea), 
                        "Shrunken_lFC" = as.data.frame(res.shr)$log2FoldChange,
                        "Shrunken_lFCSE" = as.data.frame(res.shr)$lfcSE)
  
  # 6. Gene annotation to be included 
  DEA_complete$gene_id <- rownames(DEA_complete)
  
  DEA_complete <- plyr::join(x=DEA_complete, 
             y= gene_annotations, 
             by = "gene_id", 
             type = "left", 
             match = "first")
  
  # 7. Volcano plot
  pdf(file = file.path(res.path,paste("Volcano_",comparison,suffix.fig.files,".pdf",sep="")), width=16, height=10)
  p <- EnhancedVolcano(toptable = DEA_complete,
                       lab= DEA_complete$gene_name,
                       x ='Shrunken_lFC',
                       y= 'padj',
                       xlab = "Shrunken fold change",
                       ylab = "-Log10(adj pval)",
                       title = comparison,
                       subtitle= "Differential Expression - IKba",
                       pCutoff = 0.05,
                       max.overlaps = 100,
                       labSize=3,
                       drawConnectors = TRUE,
                       FCcutoff = fc_threshold,
                       caption = paste("FC cutoff: ", fc_threshold,"; adj p-val cutoff: 0.05",sep=""),
                       legendLabels=c("Not sign.", "Shrunken FC", "adj pval","adj pval & shrunken FC"),
                       legendPosition = "right",
                       legendIconSize = 5.0,
                       legendLabSize = 10)
  print(p)
  dev.off()
  
  
  # 8. Add raw and normalized counts for those samples
  if(is.null(samples))
  {
    samples <- unlist(strsplit(x = comparison, split = "_"))[c(2,4)]
  }
  
  DEA_complete <- cbind(DEA_complete,
                        counts(DESeqSet)[,which(colData(DESeqSet)$condition %in% samples)],  # Raw
                        assay(norm.data)[,which(colData(DESeqSet)$condition %in% samples)])  # Norm
  
  # 7. Sort results by p-adj
  DEA_complete <- DEA_complete[order(DEA_complete$padj),]
  
  # 8. Store results in an excel file
  hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
  write.xlsx(DEA_complete,file = file.path(res.path,paste("DEA_All_genes_",comparison,suffix.fig.files,".xlsx",sep="")), headerStyle = hs)
  
  # 8. Return complete dataframe
  return(DEA_complete)
}



