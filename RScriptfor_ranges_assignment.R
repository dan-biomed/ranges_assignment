#lecture slides #132-139, where you have to do 2 things: 
#1. make a R markdown file for the workflow that replicates these slides. 
#2. make a new file with a new "variants" column appended to it.



source("http://bioconductor.org/biocLite.R")
biocLite()
library(GenomicRanges)
library(BiocInstaller)
biocLite("GenomicFeatures")

library(rtracklayer)

#this downloads the annotations of the mouse reference genome mm10, and puts in a transcripts database txdb.
biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene")
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
# there are many annotations accessible using accessors such as:
genes(txdb)


#imports the data file containing chromosome 1 into dbspn137
getwd()
setwd("\\Users\\dan\\Desktop\\BCB546X-Spring2017\\bds-files\\chapter-09-working-with-range-data")
dbsnp137  <- import("mm10_snp137_chr1_trunc.bed.gz")

#we need the annotation info of chromosome 1 from the reference mm10, so first collapse all overlapping exons in chr1
collapsed_exons <- reduce(exons(txdb), ignore.strand=TRUE)
chr1_collapsed_exons <- collapsed_exons[seqnames(collapsed_exons) == "chr1"]

#ready to pull out  variants that overlap exons on chromosome 1, but first, need to inspect if a variant has a width of 0,
#since if it does then we cannot find its overlap with exon ranges.
summary(width(dbsnp137))

#since the minimum shows 0, to correct this:
dbsnp137_resized <- dbsnp137
zw_i <- width(dbsnp137_resized) == 0
dbsnp137_resized[zw_i] <- resize(dbsnp137_resized[zw_i], width=1)

#now to pull out those variants that overlap exons on chromosome 1 by creating a hits object:
hits <- findOverlaps(dbsnp137_resized, chr1_collapsed_exons, 
                     ignore.strand=TRUE)
hits

#and determine the number of variants and the proportion of variants that are exonic:
length(unique(queryHits(hits)))
length(unique(queryHits(hits)))/length(dbsnp137_resized)
#results show that much more variation appear in introns

#We can also use the countOverlaps() function to find the number of variants per exon and store it it var_counts
#(note we have to reverse the order of the query since we're finding values per exon now)
var_counts <- countOverlaps(chr1_collapsed_exons, dbsnp137_resized, ignore.strand=TRUE)
var_counts

#and we can append this to our GRanges object that includes exons:
chr1_collapsed_exons$num_vars <- var_counts

#to export 
export(chr1_collapsed_exons, con="chr1_collapsed_exons.bed",
       format="bed")
chr1_collapsed_exons

