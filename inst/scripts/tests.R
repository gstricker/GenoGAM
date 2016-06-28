#################################
## A test script new functions ##
#################################

## This script should set up a working environment to test if functions
## work correctly as intended to. It's not testing for functionality!
## That's what unittests are for.

## load libraries if necessary
library(GenoGAM)

## test data
dir <- "/s/project/coreProm/data/sacCer2"
bppk <- 100
chunksize <- bppk*50
ov <- bppk*10

## read data
ggd <- GenoGAMDataSet("config.txt", chunkSize = chunksize, overhangSize = ov,
                      design = ~ s(x) + s(x, by = tfiib), directory = dir)

## or load existing data (if a fit is needed)
load("/s/project/coreProm/Michi/rdata/results_Rpb3.RData")

## test function: computeRegionSignificance
library(rtracklayer)
gffPath <- "/s/genomes/sacCer/S288C_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.pure.gff"
gff <- import(gffPath)
gff <- gff[gff$type == "gene",]
gff$name <- ifelse(!is.na(gff$gene), gff$gene, gff$Name)
seqlevels(gff)[17] <- c("chrM")
genes <- gff[, c("ID", "name")]

diffregions <- computeRegionSignificance(fit, genes)
