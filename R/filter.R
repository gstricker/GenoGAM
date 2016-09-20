#' Merge consecutive ranges in a GRanges object
#'
#' Merges consecutive range intervalls in a GRanges object to one range given a certain tolerance intervall
#'
#' @param ranges A \code{GRanges} object
#' @param overlap The overlap in basepairs. All consecutive ranges which differ by this value
#' and less will be merged. E.g. a negative number -x means an actual overlap of x+1 basepairs. A positive
#' number x means a gap of x-1 basepairs. Default is zero, which merges all consecutive ranges with at least
#' one mutual position.
#' @return A |code{GRanges} object of size <= the initial ranges object
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
mergeRanges <- function(ranges, overlap = 0) {
    mergedGRanges <- sapply(seqlevels(ranges), function(y) {
        subgr <- ranges[seqnames(ranges) == y,]
        distances <- start(subgr)[-1] - end(subgr)[-length(subgr)]
        ov <- which(distances <= overlap)
        if(length(ov) == 0) {
            return(subgr)
        }
        mergedStarts <- setdiff(ov, ov + 1)
        mergedEnds <- setdiff(ov + 1, ov)
        
        res <- c(subgr[-c(ov, ov + 1),],
                 GenomicRanges::GRanges(y, IRanges(start(subgr[mergedStarts,]),
                                                   end(subgr[mergedEnds,]))))
        return(sort(res))
    })
    return(do.call("c", unname(mergedGRanges)))
}

#' A filter function for |code{GenoGAMDataSet}
#'
#' A function to filter the |code{GenoGAMDataSet} by the sum or mean of counts
#'
#' @param ggd A \code{GenoGAMDataSet} object
#' @param threshold A value for the mean or sum of counts, which will be used to filter on basepair level.
#' By default it is taken as median + 3*MAD
#' @param windowsize The sliding window size
#' @param mode Should the sum or the mean of counts be used?
#' @return A \code{GRanges} object containing the filtered regions
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
compute_filter <- function(ggd, threshold = NULL, windowsize = 201, mode = c("sum", "mean")) {

    mode <- match.arg(mode)
    
    if(windowsize%%2 == 0) {
        windowsize <- windowsize + 1
    }
    tileSize <- tileSettings(ggd)$tileSize

    ## sliding window of sums
    sumsList <- lapply(colnames(ggd), function(y) {
        if(mode == "sum") {
            ans <- runsum(assay(ggd)[,y], windowsize, endrule = "constant")
        }
        if(mode == "mean") {
            ans <- runmean(assay(ggd)[,y], windowsize, endrule = "constant")
        }
        return(ans)
    })
    sums <- rowSums(data.frame(sumsList))

    ## compute threshold as 3 times the MAD from median when not provided
    if(is.null(threshold)) {
        sumsMedian <- median(sums)
        sumsMAD <- mad(sums)
        threshold <- sumsMedian + 3*sumsMAD
        futile.logger::flog.info(paste("Threshold estimated at", threshold))
    }

    ## find ranges complying with the threshold
    indx <- which(sums >= threshold)
    poi <- rowRanges(ggd)[indx,]
    diffs <- abs(diff(pos(poi)))
    breaks <- which(diffs > tileSize)

    starts <- c(start(poi[1,]), start(poi[breaks + 1,]))
    ends <- c(start(poi[breaks,]), start(poi[length(poi),]))
    chroms <- c(seqnames(poi[breaks,]), seqnames(poi[length(poi),]))
    gr <- GenomicRanges::GRanges(chroms, IRanges(starts, ends))

    ## correct for ranges shorter than tile size
    shortTiles <- which(width(gr) < tileSize)
    gr[shortTiles,] <- flank(gr[shortTiles,], width = ceiling(tileSize/2), both = TRUE)

    ## correct ranges outside of chromosome due to resize of short tiles
    negRanges <- which(start(gr) < 1)
    end(gr[negRanges,]) <- end(gr[negRanges,]) - start(gr[negRanges,]) + 1
    start(gr[negRanges,]) <- 1
    posRanges <- which(end(gr) > seqlengths(ggd)[as.character(seqnames(gr))])
    start(gr[posRanges,]) <- start(gr[posRanges,]) - end(gr[posRanges,]) + seqlengths(ggd)[names(posRanges)]
    end(gr[posRanges,]) <- seqlengths(ggd)[names(posRanges)]

    seqlevels(gr, force = TRUE) <- seqlevelsInUse(gr)
    res <- mergeRanges(gr, tileSize)
    return(res)
}

#' A filter function for |code{GenoGAMDataSet}
#'
#' A function to filter the |code{GenoGAMDataSet} by the sum or mean of counts to significantly reduce
#' the amount of models to compute
#'
#' @param ggd A \code{GenoGAMDataSet} object
#' @param threshold A value for the mean or sum of counts, which will be used to filter on basepair level.
#' By default it is taken as median + 3*MAD
#' @param windowsize The sliding window size. Should be an odd value.
#' @param mode Should the sum or the mean of counts be used?
#' @return A \code{GenoGAMDataSet} object containing the filtered regions
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
filterData <- function(ggd, threshold = NULL, windowsize = 201, mode = c("sum", "mean")) {

    futile.logger::flog.info("Filtering dataset for enriched regions")
    mode <- match.arg(mode)
    res <- compute_filter(ggd, threshold, windowsize, mode)
    prettyDrop <- format(nrow(ggd) - sum(width(res)), big.mark = ",", scientific = FALSE)
    prettyLeft <- format(sum(width(res)), big.mark = ",", scientific = FALSE)
    ans <- ggd[res]
    futile.logger::flog.info(paste("DONE. A total of", prettyDrop, "positions were dropped,", prettyLeft, "are left."))
    return(ans)
}
