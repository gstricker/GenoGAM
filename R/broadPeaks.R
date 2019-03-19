#' main higher function to call broad peaks. Passes arguments accordingly
#' to lower level functions
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callBroadPeaks <- function(iter, grid, fits, se, rowRanges, range, smooth, background,
                            maxgap, cutoff, is_split, is_hdf5) {
    
    if(is_split) {
        if(is_hdf5) {
            res <- .callBroadPeaks_split_hdf5(iter, grid, fits, se, rowRanges, range,
                                              smooth, background, maxgap, cutoff)
        }
        else {
            res <- .callBroadPeaks_split(iter, grid, fits, se, rowRanges, range,
                                         smooth, background, maxgap, cutoff)
        }
    }
    else {
        if(is_hdf5) {
            res <- .callBroadPeaks_hdf5(iter, grid, fits, se, rowRanges, range,
                                        smooth, background, maxgap, cutoff)
        }
        else {
            res <- .callBroadPeaks_default(iter, grid, fits, se, rowRanges, range,
                                           smooth, background, maxgap, cutoff)
        }
    }
    return(res)
}


#' Function to call broad peaks on default GenoGAM data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callBroadPeaks_default <- function(iter, grid, fits, se, rowRanges, range, smooth,
                                    background, maxgap, cutoff) {
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices 
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges, r))
    startPos <- idx[1]

    ## compute zscore for all positions
    zscore <- (fits[idx, sx] - mu0)/(sqrt(se[idx, sx]^2 + var0))
    pvals <- -pnorm(-zscore, log.p = TRUE)
    pos <- start(r):end(r)

    signifPos <- pos[pvals >= -log(cutoff)]
    if(length(signifPos) == 0) {
        res <- data.table::data.table(seqnames = character(), start = integer(),
                                      end = integer(), score = double(),
                                      meanSignal = double(), fdr = double())
    }
    else {
        diffs <- abs(diff(signifPos))
        breaks <- which(diffs > maxgap)
        
        starts <- c(signifPos[1], signifPos[breaks + 1])
        ends <- c(signifPos[breaks], signifPos[length(signifPos)])
        chrom <- GenomeInfoDb::seqnames(r)
        gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(starts, ends))
        gp <- GenomicRanges::GPos(gr)

        ## non-significant position regions  < maxgap are incorporated into the
        ## broad peak. Thus we have more positions than signifPos and have to match
        ## and normalize them to start with 1 in order to use them as an index.
        indx <- match(pos(gp), pos) - pos[1] + 1
        gp$zscore <- zscore[indx]
        gp$pval <- pvals[indx]
        gp$estimate <- exp(fits[indx, sx])
        gp$region <- S4Vectors::subjectHits(IRanges::findOverlaps(gp, gr))
        
        ## compute significance for regions just like
        ## in differential binding: 1. Hochberg for a region-wise pvalue
        ## 2. BH for FDR of the regions between each other. 
        res <- data.table::data.table(as.data.frame(gr))
        dt <- data.table::data.table(as.data.frame(gp))
        pv <- dt[, min(p.adjust(exp(-pval), method="hochberg")), by = region]
        fpb <- dt[, mean(estimate), by = region]
        res$score[pv$region] <- pv$V1
        res$meanSignal[fpb$region] <- fpb$V1
        res$fdr = p.adjust(res$score, method="BH")
        res$score <- -log(res$score)
    }

    return(res)
}

#' Function to call broad peaks on HDF5 backend data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callBroadPeaks_hdf5 <- function(iter, grid, fits, se, rowRanges, range, smooth,
                                  background, maxgap, cutoff) {
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges, r))

    ## compute zscore for all positions
    ## drop = FALSE makes sure, that both matrices are not converted to numeric
    ## vectors prior to subsetting and computation, keeping memory load low
    ## subset afterwards (pnorm can't deal with DelayedMatrix objects)
    zscore <- (fits[idx, sx, drop = FALSE] - mu0)/
        (sqrt(se[idx, sx, drop = FALSE]^2 + var0))
    pvals <- -pnorm(-zscore[,1], log.p = TRUE)
    pos <- start(r):end(r)

    signifPos <- pos[as.numeric(pvals) >= -log(cutoff)]
    if(length(signifPos) == 0) {
        res <- data.table::data.table(seqnames = character(), start = integer(),
                                      end = integer(), score = double(),
                                      meanSignal = double(), fdr = double())
    }
    else {
        diffs <- abs(diff(signifPos))
        breaks <- which(diffs > maxgap)

        starts <- c(signifPos[1], signifPos[breaks + 1])
        ends <- c(signifPos[breaks], signifPos[length(signifPos)])
        chrom <- GenomeInfoDb::seqnames(r)
        gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(starts, ends))
        gp <- GenomicRanges::GPos(gr)

        ## non-significant position regions  < maxgap are incorporated into the
        ## broad peak. Thus we have more positions than signifPos and have to match
        ## and normalize them to start with 1 in order to use them as an index.
        indx <- match(pos(gp), pos) - pos[1] + 1
        gp$zscore <- zscore[indx]
        gp$pval <- pvals[indx]
        gp$estimate <- exp(fits[indx, sx, drop = FALSE])
        gp$region <- S4Vectors::subjectHits(IRanges::findOverlaps(gp, gr))

        ## compute significance for regions just like
        ## in differential binding: 1. Hochberg for a region-wise pvalue
        ## 2. BH for FDR of the regions between each other. 
        res <- data.table::data.table(as.data.frame(gr))
        dt <- data.table::data.table(as.data.frame(gp))
        pv <- dt[, min(p.adjust(exp(-pval), method="hochberg")), by = region]
        fpb <- dt[, mean(estimate), by = region]
        res$score[pv$region] <- pv$V1
        res$meanSignal[fpb$region] <- fpb$V1
        res$fdr = p.adjust(res$score, method="BH")
        res$score <- -log(res$score)
    }
    return(res)
}

#' Function to call broad peaks on split data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callBroadPeaks_split <- function(iter, grid, fits, se, rowRanges, range, smooth,
                                  background, maxgap, cutoff) {
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices
    chr <- as.character(GenomeInfoDb::seqnames(r))
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges[[chr]], r))

    ## compute zscore for all positions
    zscore <- (fits[[chr]][idx, sx] - mu0)/(sqrt(se[[chr]][idx, sx]^2 + var0))
    pvals <- -pnorm(-zscore, log.p = TRUE)
    pos <- start(r):end(r)
    
    signifPos <- pos[pvals >= -log(cutoff)]
    if(length(signifPos) == 0) {
        res <- data.table::data.table(seqnames = character(), start = integer(),
                                      end = integer(), score = double(),
                                      meanSignal = double(), fdr = double())
    }
    else {
        diffs <- abs(diff(signifPos))
        breaks <- which(diffs > maxgap)
        
        starts <- c(signifPos[1], signifPos[breaks + 1])
        ends <- c(signifPos[breaks], signifPos[length(signifPos)])
        chrom <- GenomeInfoDb::seqnames(r)
        gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(starts, ends))
        gp <- GenomicRanges::GPos(gr)

        ## non-significant position regions  < maxgap are incorporated into the
        ## broad peak. Thus we have more positions than signifPos and have to match
        ## and normalize them to start with 1 in order to use them as an index.
        indx <- match(pos(gp), pos)
        gp$zscore <- zscore[indx]
        gp$pval <- pvals[indx]
        gp$estimate <- exp(fits[[chr]][indx, sx])
        gp$region <- S4Vectors::subjectHits(IRanges::findOverlaps(gp, gr))

        ## compute significance for regions just like
        ## in differential binding: 1. Hochberg for a region-wise pvalue
        ## 2. BH for FDR of the regions between each other. 
        res <- data.table::data.table(as.data.frame(gr))
        dt <- data.table::data.table(as.data.frame(gp))
        pv <- dt[, min(p.adjust(exp(-pval), method="hochberg")), by = region]
        fpb <- dt[, mean(estimate), by = region]
        res$score[pv$region] <- pv$V1
        res$meanSignal[fpb$region] <- fpb$V1
        res$fdr = p.adjust(res$score, method="BH")
        res$score <- -log(res$score)
    }

    return(res)
}

#' Function to call broad peaks on split HDF5 backend data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.callBroadPeaks_split_hdf5 <- function(iter, grid, fits, se, rowRanges, range, smooth,
                                  background, maxgap, cutoff) {
    r <- range[grid[iter, 1]] ## the range
    sx <- smooth[grid[iter, 2]]
    futile.logger::flog.debug(paste0("Calling peaks in region ", as.character(r)))

    mu0 <- background[sx, 'mu0']
    var0 <- background[sx, 'var0']

    ## find the region indices
    chr <- as.character(GenomeInfoDb::seqnames(r))
    idx <- S4Vectors::queryHits(IRanges::findOverlaps(rowRanges[[chr]], r))

    ## compute zscore for all positions
    ## drop = FALSE makes sure, that both matrices are not converted to numeric
    ## vectors prior to subsetting and computation, keeping memory load low
    ## subset afterwards (pnorm can't deal with DelayedMatrix objects)
    zscore <- (fits[[chr]][idx, sx, drop = FALSE] - mu0)/
        (sqrt(se[[chr]][idx, sx, drop = FALSE]^2 + var0))
    pvals <- -pnorm(-zscore[,1], log.p = TRUE)
    pos <- start(r):end(r)
    
    signifPos <- pos[as.numeric(pvals) >= -log(cutoff)]
    if(length(signifPos) == 0) {
        res <- data.table::data.table(seqnames = character(), start = integer(),
                                      end = integer(), score = double(),
                                      meanSignal = double(), fdr = double())
    }
    else {
        diffs <- abs(diff(signifPos))
        breaks <- which(diffs > maxgap)
        
        starts <- c(signifPos[1], signifPos[breaks + 1])
        ends <- c(signifPos[breaks], signifPos[length(signifPos)])
        chrom <- GenomeInfoDb::seqnames(r)
        gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(starts, ends))
        gp <- GenomicRanges::GPos(gr)

        ## non-significant position regions  < maxgap are incorporated into the
        ## broad peak. Thus we have more positions than signifPos and have to match
        ## and normalize them to start with 1 in order to use them as an index.
        indx <- match(pos(gp), pos)
        gp$zscore <- zscore[indx]
        gp$pval <- pvals[indx]
        gp$estimate <- exp(fits[[chr]][indx, sx])
        gp$region <- S4Vectors::subjectHits(IRanges::findOverlaps(gp, gr))

        ## compute significance for regions just like
        ## in differential binding: 1. Hochberg for a region-wise pvalue
        ## 2. BH for FDR of the regions between each other. 
        res <- data.table::data.table(as.data.frame(gr))
        dt <- data.table::data.table(as.data.frame(gp))
        pv <- dt[, min(p.adjust(exp(-pval), method="hochberg")), by = region]
        fpb <- dt[, mean(estimate), by = region]
        res$score[pv$region] <- pv$V1
        res$meanSignal[fpb$region] <- fpb$V1
        res$fdr = p.adjust(res$score, method="BH")
        res$score <- -log(res$score)
    }
    
    return(res)
}
