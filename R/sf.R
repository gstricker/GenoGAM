##########################
## Everything DESeq related
###############################

#' computeSizeFactors
#'
#' The function computes the size factors for given factor groups based on
#' the DESeq2 package.
#' 
#' @param ggd A GenoGAMDataSet object.
#' @param factorGroups A list of grouped IDs (same as the colnames of
#' the GenoGAMDataSet object). Each element of the list represents
#' a group of samples within which size factors are computed.
#' If NULL all samples are regarded to belong to one group.
#' Size factors are not computed between groups.
#' @return An updated GenoGAMDataSet object.
#' @examples
#' ggd <- makeTestGenoGAMDataSet()
#' ggd <- computeSizeFactors(ggd)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
computeSizeFactors <- function(ggd, factorGroups = NULL) {

    futile.logger::flog.info("Computing size factors")

    if(all(dim(ggd) == c(0, 0))) return(ggd)

    ## generate input
    sumMatrix <- sum(ggd)[[1]]
    if(is.null(factorGroups)) factorGroups <- list(colnames(ggd))
    
    ## compute sizeFactors
    sf <- NULL
        for(elem in factorGroups) {
            dds <- .normalize(sumMatrix[,elem], factor(elem))
            sf <- c(sf, log(DESeq2::sizeFactors(dds)))
        }

    ## add to GenomicTiles
    slot(ggd, "sizeFactors") <- sf

    futile.logger::flog.info("DONE")
    return(ggd)
}

#' Compute size factors.
#'
#' @param countMatrix The matrix of counts for each tile and each sample.
#' @param factors The size factors.
#' @return A DESeqDataSet object.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.normalize <- function(countMatrix, factors) {
    factors <- as.factor(factors)
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = countMatrix,
        colData = DataFrame(condition = factors),
        design = ~condition)
    dds <- DESeq2::estimateSizeFactors(dds)
    return(dds)
}

#' Perform DESeq.
#'
#' @param countMatrix The matrix of counts for each tile and each sample.
#' @param factors The size factors.
#' @return A vector of log-pvalues.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.deseq <- function(countMatrix, factors) {
    factors <- as.factor(factors)
        
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = countMatrix,
        colData = DataFrame(condition = factors),
        design = ~condition)

    dds <- DESeq2::DESeq(dds)
    ans <- DESeq2::results(dds)
    res <- ifelse(ans$log2FoldChange > 0,  log10(ans$pvalue/2), 0)
    res[is.na(res)] = 0
    return(res)
}
