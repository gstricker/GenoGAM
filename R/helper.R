.getVars <- function(formula, type = c("by", "covar")) {
    type <- match.arg(type)
    chFormula <- as.character(formula)
    variables <- gsub(" ", "", strsplit(chFormula[length(chFormula)],"\\+")[[1]])
    r <- switch(type,
                by = "by=(.+)\\)",
                covar = "s\\((.+)\\)")
    if(type == "covar") variables <- gsub(",.+", ")", variables)
    temp <- sapply(sapply(variables, function(m) regmatches(m,regexec(r,m))), function(y) y[2])
    res <- gsub(",(.+)", "", temp)
    return(res)
}

#' #' Update Formula with a specific penalization parameter lambda
#' Not used at the moment
#'
#' @param formula A formula object
#' @param lambda A vector of lambda values
#' @return Updated formula object
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.updateFormula <- function(formula, knots, m) {
    covars <- .getVars(formula, "covar")
    chFormula <- as.character(formula)
    variables <- strsplit(chFormula[2],"\\+")[[1]]
    variables[covars == "x"] <- gsub("x", "pos", variables[covars == "x"])
       
    variables <- gsub(")", paste0(', bs = "ps", k = ', knots ,', m = ', m ,')'),
                                variables)
    variablesCombined <- paste0(variables,collapse = " + ")
    reassembledFormula <- as.formula(paste("value ~ offset(offset) + ", variablesCombined))
    return(reassembledFormula)
}

#' Subset the object according to certain position in the coords.
#' 
#' @param object A data.frame alike object.
#' @param A vector of length two, specifying the start and end sites to be trimmed to
#' @return A trimmed object
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.untile <- function(object, coords = NULL){
    if(is.null(coords)) coords <- c(start(metadata(object)$chunks), end(metadata(object)$chunks))
    start <- min(object$pos)
    end <- max(object$pos)
    overhang <- c(which(object$pos < coords[1]), which(object$pos >= coords[2]))
    return(object[-overhang,])
}

#' Create a melted version of GenomicTiles,adding the experimentMatrix columns
#' and melting by samplenames.
#'
#' @param gtiles, A DataFrameList alike object.
#' @param settings A list of tile settings.
#' @param colData The sample specific values, to be added to the melted data.frame.
#' The rownames of colData should comply with the sample names in the config data.frame.
#' @param experimentMatrix A matrix object representing the experiment design.
#' @return A melted data.frame with additional columns.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.meltGTile <- function(gtiles, experimentDesign, sf, formula) {
    res <- lapply(gtiles, function(gt) {

        ids <- rownames(experimentDesign)
        ## melt
        melted <- reshape2::melt(as.data.frame(gt), measure.vars = ids,
             variable.name = "ID")

        ## add offset
        if(!("offset" %in% names(melted))) {
            melted$offset <- 0
        }

        rle <- Rle(melted$ID)
        sf <- sf[match(runValue(rle), names(sf))] ## put them in right order
        melted$offset <- melted$offset + rep(sf, runLength(rle))
        
        ## add experimentDesign columns
        designVars <- .getVars(formula, "by")
        designCols <- experimentDesign[names(sf),na.omit(designVars), drop = FALSE]
        for(col in names(designCols)) {
            melted[[col]] <- rep(designCols[[col]], runLength(rle))
        }
        return(melted)
    })
    names(res) <- names(gtiles)
    return(res)
}

#' Merge consecutive ranges in a GRanges object
#'
#' Merges consecutive range intervalls in a GRanges object to one range given a certain tolerance intervall
#'
#' @param ranges A |code{GRanges} object
#' @param overlap The overlap in basepairs. All consecutive ranges which differ by this value
#' and less will be merged. E.g. a negative number -x means an actual overlap of x+1 basepairs. A positive
#' number x means a gap of x-1 basepairs. Default is zero, which merges all consecutive ranges with at least
#' one mutual position.
#' @return A |code{GRanges} object of size <= the initial ranges object
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
mergeRanges <- function(ranges, overlap = 0) {
    mergedGRanges <- sapply(seqlevels(ranges), function(y) {
        subgr <- ranges[seqnames(ranges) == y,]
        distances <- start(subgr)[-1] - end(subgr)[-length(subgr)]
        ov <- which(distances <= overlap)
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
#' @param ggd A |code{GenoGAMDataSet} object
#' @param threshold A value for the mean or sum of counts, which will be used to filter on basepair level.
#' By default it is taken as median + 3*MAD
#' @param windowsize The sliding window size
#' @param mode Should the sum or the mean of counts be used?
#' @return A |code{GRanges} object containing the filtered regions
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
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

    res <- mergeRanges(gr, tileSize)
    return(res)
}

filter <- function(ggd, threshold = NULL, windowsize = 201, mode = c("sum", "mean")) {
   
    mode <- match.arg(mode)
    res <- compute_filter(ggd, threshold, windowsize, mode)
    return(ggd[res])
}
    
