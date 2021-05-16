setwd("/media/ggj/Files/scATAC/snapATAC/COL10/")

library(ATACseqQC)
## input the bamFile from the ATACseqQC package 
bamfile <- "tmp/human.bam"
bamfile.labels <- gsub(".bam", "", basename(bamfile))

#Estimate the library complexity
bamQC(bamfile, outPath="./")
estimateLibComplexity(readsDupFreq(bamfile))

## generate fragement size distribution
fragSize <- fragSizeDist(bamfile, bamfile.labels)
