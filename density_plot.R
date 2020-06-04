
#source("https://bioconductor.org/biocLite.R")
biocLite("karyoploteR")
library(karyoploteR)

custom.genome <- toGRanges(data.frame(chr=c("NC_035780.1","NC_035781.1","NC_035782.1","NC_035783.1","NC_035784.1","NC_035785.1",
                                            "NC_035786.1", "NC_035787.1","NC_035788.1","NC_035789.1"), 
                                      start=c(1,1,1,1,1,1,1,1,1,1), 
                                      end=c(65668440,61752955,77061148,59691872,98698416,51258098,57830854,75944018,104168038,32650045)))
custom.genome_b <- toGRanges(data.frame(chr=c("NC_035780.1"), 
                                        start=c(1), 
                                        end=c(10000)))

kp <- plotKaryotype(genome = custom.genome)
kp_b <- plotKaryotype(genome = custom.genome_b)

#read the data with regioneR's toGRanges
cnv <- toGRanges("/Users/tejashree/Documents/Projects/cnv/delly/oysterduplicate_sort.bed")
cnv_fil <- toGRanges("/Users/tejashree/Documents/Projects/cnv/scripts/cnv/cvir_filtered_dups.bed")
exon <- toGRanges("/Users/tejashree/Documents/Projects/cnv/annot/exon.bed")

#ifi44 dup distribution against ifi44 exon distribution
ifi44_dist <- toGRanges("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/ifi44_dup.bed")
ifi44_exons <- toGRanges("/Users/tejashree/Documents/Projects/cnv/annot/ifi44_exon.bed")

#gimap dup distribution against ifi44 exon distribution
gimap_dist <- toGRanges("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/gimap_dup.bed")
gimap_exons <- toGRanges("/Users/tejashree/Documents/Projects/cnv/annot/gimap_exon.bed")

#kpPlotDensity(kp, data=cnv, window.size = 0.5e6, col = "blue")
kpPlotDensity(kp, data=cnv, window.size = 300000, col = "red", r0=0.5, r1=1)
kpPlotDensity(kp, data=cnv_fil, window.size = 300000, col = "red", r0=0.5, r1=1)
kpPlotDensity(kp, data=exon, window.size = 300000, col = "blue", r0=0.5, r1=0)

kpPlotDensity(kp_b, data=cnv_fil, window.size = 100, col = "red", r0=0.5, r1=1)
kpPlotDensity(kp_b, data=exon, window.size = 100, col = "blue", r0=0.5, r1=0)


#Map IFI44 dups and ifi44 exons on one plot
#IFI44 dups
kpPlotDensity(kp, data=ifi44_dist, window.size = 500, col = "green", r0=0.5, r1=1)
kpPlotDensity(kp, data=ifi44_exons, window.size = 500, col = "orchid", r0=0.5, r1=0)

#Map gimap dups and gimap exons on one plot
#GIMAP dups
kpPlotDensity(kp, data=gimap_dist, window.size = 500, col = "green", r0=0.5, r1=1)
kpPlotDensity(kp, data=gimap_exons, window.size = 500, col = "orchid", r0=0.5, r1=0)

