###############################
## size factors correction
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
#' @return A GenoGAMDataSet object, where the sizeFactors slot is updated.
#' @examples
#' ggd <- makeTestGenoGAMDataSet()
#' ggd <- computeSizeFactors(ggd)
#' sizeFactors(ggd)
#' groups <- list(c("wt_1", "wt_2"), c("mutant_1", "mutant_2"))
#' ggd <- computeSizeFactors(ggd, factorGroups = groups)
#' sizeFactors(ggd)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
computeSizeFactors <- function(ggd, factorGroups = NULL) {

    futile.logger::flog.info("Computing size factors")

    if(all(dim(ggd) == c(0, 0))) return(ggd)
    if(is.null(colnames(ggd))) {
        futile.logger::flog.error("GenoGAMDataSet doesn't have column names. No size factors computed.")
        return(ggd)
    }

    ## generate input
    sumMatrix <- getCountMatrix(ggd)
    if(is.null(factorGroups)) {
        factorGroups <- list(colnames(ggd))
    }

    futile.logger::flog.debug("Using the following factor groups:", factorGroups, capture = TRUE)
    
    ## compute sizeFactors
    sf <- NULL
    for(elem in factorGroups) {
        dds <- .normalize(sumMatrix[,elem, drop = FALSE], factor(elem))
        sf <- c(sf, log(DESeq2::sizeFactors(dds)))
    }
    idx <- which(!(colnames(ggd) %in% names(sf)))
    if(length(idx) > 0) {
        zeros <- rep(0, length(idx))
        names(zeros) <- colnames(ggd)[idx]
        sf <- c(sf, zeros)
        sf <- sf[match(colnames(ggd), names(sf))]
        names(sf) <- colnames(sf)
    }
    
    ## add to GenomicTiles
    futile.logger::flog.debug("The log-size-factors were computed as:", sf, capture = TRUE)
    slot(ggd, "sizeFactors") <- sf

    futile.logger::flog.info("DONE")
    return(ggd)
}

#' Compute size factors.
#'
#' @param countMatrix The matrix of counts for each tile and each sample.
#' @param factors The size factors.
#' @return A DESeqDataSet object.
#' @noRd
.normalize <- function(countMatrix, factors) {
    factors <- as.factor(factors)
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = countMatrix,
        colData = S4Vectors::DataFrame(condition = factors),
        design = ~condition)
    dds <- DESeq2::estimateSizeFactors(dds)
    return(dds)
}

#' Perform DESeq.
#'
#' @param countMatrix The matrix of counts for each tile and each sample.
#' @param factors The size factors.
#' @return A vector of log-pvalues.
#' @noRd
.deseq <- function(countMatrix, factors) {
    if(nrow(factors) == (ncol(factors) + 1)) {
        futile.logger::flog.info("The design matrix has the same number of samples and coefficients. The top regions for cross-validation will be determined by the highest sum of counts.")
        ## The -1 is necessary to be consistent with the ordering of p-values calculated in
        ## the 'else' clause. P-values go from smallest (most significant) to biggest
        ## (least significant). 
        res <- -1 * rowSums(countMatrix)
    }
    else {
        form <- as.formula(paste0("~", paste(names(factors), collapse = "+")))
        
        dds <- DESeq2::DESeqDataSetFromMatrix(
                           countData = countMatrix,
                           colData = factors,
                           design = form)

        dds <- DESeq2::DESeq(dds)
        ans <- DESeq2::results(dds)
        res <- ifelse(ans$log2FoldChange > 0,  log10(ans$pvalue/2), 0)
        res[is.na(res)] = 0
    }
    return(res)
}
