################################################################################
##        Title: IkBa HSC dormancy - Thambyrajah R. et al. 2024
##                "IκBα controls dormancy in Hematopoietic Stem Cells via 
##                retinoic acid during embryonic development"
## Description :  This script is to prepare a list of regions (exclusive for WT or KO)
##               for conducting motif analysis i.e. MEME suite (MEME-ChIP) recommends sequences of 500 bp length
################################################################################

require(GenomicRanges)

# Exclusive regions are listed in an excel file: This can be taken from Suppl Data T4

exclusive.list <- list("exclusive.WT" = read.xlsx(xlsxFile = "Supp_Data_T4.xlsx", sheet="Exclusive.WT.FL"),
                       "exclusive.KO" = read.xlsx(xlsxFile = "Supp_Data_T4.xlsx", sheet="Exclusive.KO.FL"))
scenario <- "FL.H3K27me3.Top0.5" # This is required in case there are different scenarios to be considered

#######################
# Create GRanges object (load functio at the end of this script)
gr.exclusive.list <- lapply(exclusive.list, function(df)  create.GRanges(df.regions = df))

# Split regions into bins of 500bp with an sliding window
bin = 500 
gr.exclusive.list.split <- lapply(gr.exclusive.list, function(l) 
  {
  res <- as.data.frame(slidingWindows(x=l, width = bin, step = bin))
  
  for(i in 1:nrow(res))
  {
    if(res[i,"width"] < 500)
      res[i,"end"] <- res[i,"end"] + (500-res[i,"width"])
  }
  return(res)
  })

# Create bed file from this
lapply(names(gr.exclusive.list.split), function(case) 
{
  sp <- gr.exclusive.list.split[[case]]  
  
  write.table(x = data.frame("chrom" = sp$seqnames,
                             "start" = sp$start,
                             "down" = sp$end),
            file=paste("Ikba_HSC_ROSHANA/BED_Exclusive/",case,"_MEMEbed_", scenario, ".bed", sep=""),
            row.names = FALSE, col.names = FALSE,
            quote = FALSE,sep="\t")
  })

# FUNCTIONS
############

create.GRanges <- function(df.regions)
{
   GRanges(seqnames = df.regions$seqnames,
            IRanges(start=df.regions$start,
                    end=df.regions$end,
                    names=df.regions$Peak_Name),
            strand=df.regions$strand)
}
