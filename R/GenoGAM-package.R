#' GenoGAM: A package providing a framework to analyse ChIP-Seq data
#' 
#' @name GenoGAM
#' @import methods
#' @import HDF5Array
#' @import rhdf5
#' @import BiocParallel
#' @import IRanges
#' @import GenomicRanges
#' @importFrom futile.logger flog.info
#' @importFrom futile.logger flog.warn
#' @importFrom futile.logger flog.error
#' @importFrom futile.logger flog.trace
#' @importFrom futile.logger flog.debug
#' @importFrom futile.logger flog.threshold
#' @importFrom stats runif rnbinom as.formula dnbinom optim na.omit pnorm p.adjust stepfun
#' @importFrom SummarizedExperiment assay assays assays<- colData colData<- SummarizedExperiment rowRanges rowRanges<-
#' @importFrom S4Vectors metadata metadata<- DataFrame Rle queryHits
#' @importFrom DelayedArray rowRanges
#' @importFrom DESeq2 design design<- sizeFactors
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqlengths
#' @importFrom GenomicAlignments seqlevelsInUse
#' @useDynLib GenoGAM
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
