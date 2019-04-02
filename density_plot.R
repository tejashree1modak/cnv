
#source("https://bioconductor.org/biocLite.R")
biocLite("karyoploteR")
library(karyoploteR)

custom.genome <- toGRanges(data.frame(chr=c("NC_035780.1","NC_035781.1","NC_035782.1","NC_035783.1","NC_035784.1","NC_035785.1",
                                            "NC_035786.1", "NC_035787.1","NC_035788.1","NC_035789.1"), 
                                      start=c(1,1,1,1,1,1,1,1,1,1), 
                                      end=c(65668440,61752955,77061148,59691872,98698416,51258098,57830854,75944018,104168038,32650045)))
kp <- plotKaryotype(genome = custom.genome)

#read the data with regioneR's toGRanges
cnv <- toGRanges("/Users/tejashree/Documents/Projects/cnv/delly/oysterduplicate_sort.bed")
exon <- toGRanges("/Users/tejashree/Documents/Projects/cnv/annot/exon.bed")

#kpPlotDensity(kp, data=cnv, window.size = 0.5e6, col = "blue")
kpPlotDensity(kp, data=cnv, window.size = 300000, col = "red", r0=0.5, r1=1)
kpPlotDensity(kp, data=exon, window.size = 300000, col = "blue", r0=0.5, r1=0)
