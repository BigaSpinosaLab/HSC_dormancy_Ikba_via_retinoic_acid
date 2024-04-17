################################################################################
##        Title   : IkBa HSC dormancy - Thambyrajah R. et al. 2024
##                "IκBα controls dormancy in Hematopoietic Stem Cells via 
##                retinoic acid during embryonic development"
##  Description : This script is for creating: annotations of called peaks per sample,
##                consensus peaksets and exclusive peaksets for WT and KO regarding
##                CUT&Tag H3K27me3 data
##        Author: María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
################################################################################

# NOTES:

# a) Peaks were called with SEACR without IgG and considering Top 0.05%. BED 
#   files with peaks are being blacklisted-filtered. THESE BED FILES ARE AVAILABLE
#   IN GEO accession no. GSE188524 and are required for this analysis

# b) Consensus peakset: peaks that are detected in at least two out of the three
#    replicates - overlapping with a max gap of 100bp

# REMARK: For peaks annotation, a TxDb object is created from gtf file (Ensembl,
#         release102 mm10)
 
# REMARK: The execution of this script results in information gathered in Suppl
#         Data T4

################################################################################
## 1. Load packages
################################################################################

require(ChIPpeakAnno)
library(rtracklayer)
require(openxlsx)

# Required for annotation
require(ChIPseeker)  
require(GenomicFeatures)  # For creating a TxDb object from gtf file
require(org.Mm.eg.db)

################################################################################
## 2. Import files with called peaks
################################################################################

# Path to where all BED files are stored (i.e. downloaded from GEO)
path.SEACR <- "../BED_files"
scenarios <- c("Top0.5") # This is required in case different scenarios are explored

#====== H3K27me3 LTHSC Fetal liver WT or KO condition ============#

samples.wt <- c("Ikba_WT_LT-HSC_1","Ikba_WT_LT-HSC_2","Ikba_WT_LT-HSC_3")
samples.ko <- c("Ikba_KO_LT-HSC_1","Ikba_KO_LT-HSC_2", "Ikba_KO_LT-HSC_3")


# Define the file pattern (SEACR criteria)
pattern <- c(".woDups.top0.5.stringent.filtered_blacklisted.bed")

# Define the function for reading the beds and preparing colnames
read.beds <- function(sample, path, pattern)
{
  a <- read.table(file = paste0(path.SEACR, sample, pattern), header=FALSE)
  colnames(a) <- c("chrom", "start", "end", "total.signal", "max.signal", "opt.region")
  return(a)
}

# Read bed files for H3K27me3 Fetal Liver: read WT and KO separated
fl.wt <- lapply(pattern, function(case){
   res <- lapply(samples.wt, function(sample) read.beds(sample = sample,path = path.SEACR, pattern = case))
   names(res) <- paste0(samples.wt,"_H3K27m3")
   return(res)
})

fl.ko <- lapply(pattern, function(case){
  res <- lapply(samples.ko, function(sample) read.beds(sample = sample,path = path.SEACR, pattern = case))
  names(res) <- paste0(samples.ko,"_H3K27m3")
  return(res)
})


names(fl.wt) <- names(fl.ko) <- scenarios

################################################################################
## 2. Make GRanges objects with all intervals
################################################################################

# Function to create GRanges objects
create.GRanges <- function(list.bedfiles)
{
  gr <- lapply(names(list.bedfiles), function(mark) 
  {
    GRanges(seqnames = list.bedfiles[[mark]]$chrom,
            IRanges(start=list.bedfiles[[mark]]$start,
                    end=list.bedfiles[[mark]]$end,
                    names=paste(mark,seq(1,nrow(list.bedfiles[[mark]])),sep="_")),
            strand="*")
  })
  names(gr) <- names(list.bedfiles)
  
  return(gr)
}

gr.fl.ko <- lapply(fl.ko, function(case) create.GRanges(list.bedfiles = case))
gr.fl.wt <- lapply(fl.wt, function(case) create.GRanges(list.bedfiles = case))

################################################################################
## 3. Find overlapping peaks (100bp of maximum gap). These will be merged
################################################################################

# Overlapping peaks

overlapping.ko <- lapply(gr.fl.ko, function(case) 
  findOverlapsOfPeaks(case[[1]],case[[2]], case[[3]],
                      ignore.strand = TRUE, 
                      connectedPeaks = "keepAll",
                      maxgap = 100)  # Maxgap value between peaks: 100bp between peaks
)

overlapping.wt <- lapply(gr.fl.wt, function(case) 
  findOverlapsOfPeaks(case[[1]],case[[2]], case[[3]],                  
    ignore.strand = TRUE, 
                      connectedPeaks = "keepAll",
                      maxgap = 100)  # Maxgap value between peaks: 100bp between peaks
)

################################################################################
## 5. Create txdb for Annotate the consensus PeakSet: just execute it once
################################################################################

# IMPORTANT: Download GTF file from Ensembl Release 102 mm10

# Create the txdb according to the same used during alignment
metadata <- data.frame(name="Resource URL",
                       value=paste0("ftp://ftp.ensemblgenomes.org/pub/","release-102/gtf/mus_musculus/"))

txdb <- makeTxDbFromGFF(file = "/Volumes/projectscomput/cancer/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.102.gtf",
                        format="gtf",
                        dataSource="Ensembl_FTP_repository",
                        organism="Mus Musculus",
                        taxonomyId=10090,
                        circ_seqs=NULL,
                        chrominfo=NULL,
                        miRBaseBuild=NA,
                        metadata=metadata,
                        dbxrefTag="gene_id")


################################################################################
## 6. Annotate the consensus PeakSet and the individual peaks
################################################################################

# Function to conduct annotation to peaks
annotating_peaks <- function(peaks)
{
  annotatePeak(peak = peaks, 
               tssRegion = c(-5000,100), 
               TxDb = txdb,
               genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                             "Downstream", "Intergenic"),  # Default genomic Annotation
               annoDb = "org.Mm.eg.db")
}

# Annotate Consensus sets for KO and WT
Consensus.ko <- lapply(overlapping.ko, function(case) annotating_peaks(peaks = case$mergedPeaks))
Consensus.wt <- lapply(overlapping.wt, function(case) annotating_peaks(peaks = case$mergedPeaks))

# Add annotations to the individual samples
Individual.samples.ko <- lapply(gr.fl.ko, function(case) {
  lapply(case, function(set) {
    p <- annotating_peaks(peaks = set)
    return(as.data.frame(p)) })
}) 

Individual.samples.wt <- lapply(gr.fl.wt, function(case) {
  lapply(case, function(set) {
    p <- annotating_peaks(peaks = set)
    return(as.data.frame(p)) })
}) 

final.results <- list("Top0.5" = c(Individual.samples.wt$Top0.5, 
                                  "Consensus.WT" = list(as.data.frame(Consensus.wt$Top0.5)),
                                 Individual.samples.ko$Top0.5, 
                                 "Consensus.KO" = list(as.data.frame(Consensus.ko$Top0.5))))

################################################################################
## 7. Find overlapping peaks (by at least 100bp) and genes between both conditions
################################################################################

# Function to create a GR from the consensus
create.GRanges.consensus <- function(consensus)
{
  gr.res <- GRanges(seqnames = consensus$seqnames,
            IRanges(start=consensus$start,
                    end=consensus$end),
            strand="*")
  return(gr.res)
}


# Overlapping peaks
overlapping.conditions.consensus <- lapply(final.results, function(case) 
  findOverlapsOfPeaks(create.GRanges.consensus(case$Consensus.WT), 
                      create.GRanges.consensus(case$Consensus.KO),
                      ignore.strand = TRUE, 
                      connectedPeaks = "keepAll",
                      maxgap = 100))  # Maxgap default value: 100bp between peaks


# Store the unique peaks (not overlapping) for further study (i.e. Motif analysis)

fl.exclusive.regions <- as.data.frame(overlapping.conditions.consensus$Top0.5$uniquePeaks)
fl.exclusive.regions$Peak_Name <- names(overlapping.conditions.consensus$Top0.5$uniquePeaks)

saveRDS(object = fl.exclusive.regions, file="../FL_H3K27me3_exclusive_regions.RDS")

################################################################################
## 8. (OPTIONAL) Venn diagram on annotated genes to consensus peaks
################################################################################

# Venn Diagram with annotated genes from corresponding consensus sets: Fetal Liver

plots <- lapply(scenarios, function(case){
  
   ggvenn::ggvenn(list("WT" = final.results[[case]]$Consensus.WT$ENTREZID, 
                           "KO"= final.results[[case]]$Consensus.KO$ENTREZID),
                      stroke_color= "white", stroke_size=0.2,
                      text_size = 3.5) +
  ggtitle(paste0("Peak calling criteria: ",case))
})


pdf(file = "../FL_H3K27me3_VennDiagram_WT_KO_annotated_genes.pdf", 
     width=8, height=6)
print(plots)
dev.off()

################################################################################
## 9. Store results
################################################################################


final.results.complete  <- lapply(scenarios, function(case) {
  
  common.fl <- as.data.frame(overlapping.conditions.consensus[[case]]$mergedPeaks)
  res <- c(final.results[[case]],
    "Common.Regions.KO.WT.FL" = list(common.fl),
    "Exclusive.WT.FL" = list(clusterProfiler::bitr(geneID = setdiff(final.results[[case]]$Consensus.WT$ENTREZID, 
                                                               final.results[[case]]$Consensus.KO$ENTREZID), 
                                              fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db, drop=FALSE)),
    "Exclusive.KO.FL" = list(clusterProfiler::bitr(geneID = setdiff(final.results[[case]]$Consensus.KO$ENTREZID, 
                                                               final.results[[case]]$Consensus.WT$ENTREZID), 
                                              fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db, drop=FALSE)))
  return(res)
})
names(final.results.complete) <- scenarios


# This information is stored in Suppl Data T4
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(x= final.results.complete$Top0.5, 
           file = "../H3K27me3_TOP0.5_criteria_PEAKS_FINAL.xlsx",
           headerStyle=hs)


