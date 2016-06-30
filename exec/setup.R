## This script should set up a working environment to test if functions
## work correctly as intended to. It's not testing for functionality!
## That's what unittests are for.

## setup environment
## =================

## load libraries if necessary
library(GenoGAM)

## test data
dir <- "/s/project/coreProm/data/sacCer2"
bppk <- 100
chunksize <- bppk*50
ov <- bppk*10

## read data
if(functionName == "qualityCheck" | functionName == "makeTiles" | functionName == "filter") {
    ggd <- GenoGAMDataSet("config.txt", chunkSize = chunksize, overhangSize = ov,
                          design = ~ s(x) + s(x, by = tfiib), directory = dir)
}

## or load existing data (if a fit is needed)
if(functionName == "computeRegionSignificance") {
    load("/s/project/coreProm/Michi/rdata/results_Rpb3.RData")
}
