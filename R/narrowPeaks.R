#' main higher function to call narrow peaks. Passes arguments accordingly
#' to lower level functions
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callNarrowPeaks <- function(iter, grid, fits, se, rowRanges, range, smooth,
                             background, is_split, is_hdf5) {
    
    if(is_split){
        if(is_hdf5) {
            res <- .callNarrowPeaks_split_hdf5(iter, grid, fits, se, rowRanges, range,
                                               smooth, background)
        }
        else {
            res <- .callNarrowPeaks_split(iter, grid, fits, se, rowRanges, range,
                                          smooth, background)
        }
    }
    else {
        if(is_hdf5) {
            res <- .callNarrowPeaks_hdf5(iter, grid, fits, se, rowRanges, range,
                                         smooth, background)
        }
        else {
            res <- .callNarrowPeaks_default(iter, grid, fits, se, rowRanges, range,
                                            smooth, background)
        }
    }
    return(res)
}

#' Function to call narrow peaks on default GenoGAM data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callNarrowPeaks_default <- function(iter, grid, fits, se, rowRanges, range,
                                     smooth, background) {
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices 
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges, r))
    startPos <- idx[1]

    ## find peaks
    region <- fits[idx, sx]
    peaks <- .findPeaks(region) + startPos - 1

    if(length(peaks) == 0) {
        dfpeaks <- data.table::data.table(seqnames = character(), pos = integer(),
                                          zscore = double(), score = double(),
                                          fdr = double(), summit = double())
    }
    else {
        zscore <- (fits[peaks, sx] - mu0)/(sqrt(se[peaks, sx]^2 + var0))
        pvals <- -pnorm(-zscore, log.p = TRUE)
        dfpeaks <- data.table::data.table(seqnames = as.character(seqnames(r)),
                                          pos = peaks,
                                          zscore = zscore, score = pvals)
    
        ## for valleys
        valleys <- .findValleys(region) + startPos - 1
        vzscore <- (fits[valleys, sx] - mu0)/(sqrt(se[valleys, sx]^2 + var0))
        vzscore <- -vzscore
        dfvalleys <- data.table::data.table(pos = valleys, zscore = vzscore)
        
        ## compute FDR
        dfpeaks <- dfpeaks[order(zscore, decreasing = TRUE),]
        dfvalleys <- dfvalleys[order(zscore, decreasing = TRUE)]
        fdr <- sapply(dfpeaks$zscore, function(y) {
            sum(dfvalleys$zscore >= y)/sum(dfpeaks$zscore >= y)
        })
        ## adjust fdr for odd length of both vectors -->
        ## Number of peaks and valleys might differ by one, causing FDR become
        ## slightly greater than 1
        fdr <- fdr*(nrow(dfpeaks)/nrow(dfvalleys))
        dfpeaks$fdr <- fdr
        dfpeaks$summit <- exp(fits[peaks, sx])
    }

    return(dfpeaks)
}

#' Function to call narrow peaks on HDF5 backend data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callNarrowPeaks_hdf5 <- function(iter, grid, fits, se, rowRanges, range,
                                  smooth, background) {
    
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges, r))
    startPos <- idx[1]

    ## find peaks
    region <- fits[idx, sx, drop = FALSE]
    peaks <- .findPeaks_hdf5(region) + startPos - 1
    
    if(length(peaks) == 0) {
        dfpeaks <- data.table::data.table(seqnames = character(), pos = integer(),
                                          zscore = double(), score = double(),
                                          fdr = double(), summit = double())
    }
    else {
        zscore <- (fits[peaks, sx, drop = FALSE] - mu0)/
            (sqrt(se[peaks, sx, drop = FALSE]^2 + var0))
        pvals <- -pnorm(-zscore, log.p = TRUE)
        dfpeaks <- data.table::data.table(chromosome = as.character(seqnames(r)),
                                          pos = peaks,
                                          zscore = as.numeric(zscore),
                                          score = as.numeric(pvals))

        ## for valleys
        valleys <- .findValleys_hdf5(region) + startPos - 1
        vzscore <- (fits[valleys, sx, drop = FALSE] - mu0)/
            (sqrt(se[valleys, sx, drop = FALSE]^2 + var0))
        vzscore <- -vzscore
        dfvalleys <- data.table::data.table(pos = valleys, zscore = as.numeric(vzscore))

        ## compute FDR
        dfpeaks <- dfpeaks[order(zscore, decreasing = TRUE),]
        dfvalleys <- dfvalleys[order(zscore, decreasing = TRUE)]
        fdr <- sapply(dfpeaks$zscore, function(y) {
            sum(dfvalleys$zscore >= y)/sum(dfpeaks$zscore >= y)
        })
        ## adjust fdr for odd length of both vectors -->
        ## Number of peaks and valleys might differ by one, causing FDR become
        ## slightly greater than 1
        fdr <- fdr*(nrow(dfpeaks)/nrow(dfvalleys))
        dfpeaks$fdr <- fdr
        dfpeaks$summit <- as.numeric(exp(fits[peaks, sx, drop = FALSE]))
    }

    return(dfpeaks)
}

#' Function to call narrow peaks on split data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callNarrowPeaks_split <- function(iter, grid, fits, se, rowRanges, range,
                                   smooth, background) {
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices
    chr <- as.character(seqnames(r))
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges[[chr]], r))

    ## find peaks
    region <- fits[[chr]][idx, sx]
    peaks <- .findPeaks(region)

    if(length(peaks) == 0) {
        dfpeaks <- data.table::data.table(seqnames = character(), pos = integer(),
                                          zscore = double(), score = double(),
                                          fdr = double(), summit = double())
    }
    else {
        zscore <- (fits[[chr]][peaks, sx] - mu0)/(sqrt(se[[chr]][peaks, sx]^2 + var0))
        pvals <- -pnorm(-zscore, log.p = TRUE)
        dfpeaks <- data.table::data.table(chromosome = chr,
                                          pos = peaks,
                                          zscore = zscore, score = pvals)

        ## for valleys
        valleys <- .findValleys(region)
        vzscore <- (fits[[chr]][valleys, sx] - mu0)/(sqrt(se[[chr]][valleys, sx]^2 + var0))
        vzscore <- -vzscore
        dfvalleys <- data.table::data.table(pos = valleys, zscore = vzscore)

        ## compute FDR
        dfpeaks <- dfpeaks[order(zscore, decreasing = TRUE),]
        dfvalleys <- dfvalleys[order(zscore, decreasing = TRUE)]
        fdr <- sapply(dfpeaks$zscore, function(y) {
            sum(dfvalleys$zscore >= y)/sum(dfpeaks$zscore >= y)
        })
        ## adjust fdr for odd length of both vectors -->
        ## Number of peaks and valleys might differ by one, causing FDR become
        ## slightly greater than 1
        fdr <- fdr*(nrow(dfpeaks)/nrow(dfvalleys))
        dfpeaks$fdr <- fdr
        dfpeaks$summit <- exp(fits[[chr]][peaks, sx])
    }

    return(dfpeaks)
}

#' Function to call narrow peaks on split HDF5 backend data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callNarrowPeaks_split_hdf5 <- function(iter, grid, fits, se, rowRanges, range,
                                        smooth, background) {
    
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices
    chr <- as.character(seqnames(r))
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges[[chr]], r))

    ## find peaks
    region <- fits[[chr]][idx, sx, drop = FALSE]
    peaks <- .findPeaks_hdf5(region)
    
    if(length(peaks) == 0) {
        dfpeaks <- data.table::data.table(seqnames = character(), pos = integer(),
                                          zscore = double(), score = double(),
                                          fdr = double(), summit = double())
    }
    else {
        zscore <- (fits[[chr]][peaks, sx] - mu0)/
            (sqrt(se[[chr]][peaks, sx]^2 + var0))
        pvals <- -pnorm(-zscore, log.p = TRUE)
        dfpeaks <- data.table::data.table(chromosome = chr,
                                          pos = peaks,
                                          zscore = as.numeric(zscore),
                                          score = as.numeric(pvals))

        ## for valleys
        valleys <- .findValleys_hdf5(region)
        vzscore <- (fits[[chr]][valleys, sx] - mu0)/
            (sqrt(se[[chr]][valleys, sx]^2 + var0))
        vzscore <- -vzscore
        dfvalleys <- data.table::data.table(pos = valleys, zscore = as.numeric(vzscore))

        ## compute FDR
        dfpeaks <- dfpeaks[order(zscore, decreasing = TRUE),]
        dfvalleys <- dfvalleys[order(zscore, decreasing = TRUE),]
        fdr <- sapply(dfpeaks$zscore, function(y) {
            sum(dfvalleys$zscore >= y)/sum(dfpeaks$zscore >= y)
        })
        ## adjust fdr for odd length of both vectors -->
        ## Number of peaks and valleys might differ by one, causing FDR become
        ## slightly greater than 1
        fdr <- fdr*(nrow(dfpeaks)/nrow(dfvalleys))
        dfpeaks$fdr <- fdr
        dfpeaks$summit <- as.numeric(exp(fits[[chr]][peaks, sx]))
    }

    return(dfpeaks)
}
