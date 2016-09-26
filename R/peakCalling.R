#' Call peaks on a GenoGAM object
#'
#' Call narrow or broad peaks on the GenoGAM fit and computing significance, respectively 
#'
#' @param fit A \code{GenoGAM} object
#' @param smooth The name of the smooth, i.e. the 'by' variables in the GenoGAMDataSet design. By default
#' the last one will be taken.
#' @param range A \code{GRanges} object specifying a range. Bt default the complete fit is taken.
#' @param peakType The type of the peak (narrow or broad). Default is narrow, see details.
#' @param threshold The significance threshold. Keep in mind that the treshold depends on the thresholdType.
#' By default this is 0.05 for 'pvalue' and 0.1 for 'fdr'.
#' @param thresholdType The threshold type. Either 'fdr'(default) or 'pvalue'. If the threshold is not provided it,
#' will be set accordingly to the thresholdType.
#' @param maxgap For broad peaks only. The maximum gap between two broad peaks, that can be tolerated in order to
#' identify both as part of one broad peak. All broad peaks with distances smaller or equal to the maxgap will be
#' merged.
#' @param cutoff A seperate threshold for broad peaks. Since pointwise pvalues are available, this threshold is used
#' to identify all significantly high positions, which then make up a broad peak.
#' @param minregion For broad peaks only. The minimum length of a broad peak. By default 1, thus catching also narrow peaks.
#' @return A data.table of identified peaks. The different columns loosely resemble the narrow and broad peak format
#' (with differentcolumn order), such that it is easy to write them to a 'narrowPeak', 'broadPeak' file. See details for
#' column description.
#' @details Note, that broad peaks don't provide a specific highest location, but a region. Whereas narrow peaks
#' provide both. However, the borders of narrow peaks are not necessarily informative. Additionally narrow peaks provide
#' a 95\% confidence interval for the position, namely 'start' and 'end', which gives a more informative uncertainty
#' measure to the peak position. Also narrow peaks provide an occupancy estimate at the peak position, while broad peaks
#' give the average occupancy accross the region.
#' The columns returned are:
#' 
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
callPeaks <- function(fit, smooth = NULL, range = NULL, peakType = c("narrow", "broad"), threshold = NULL, thresholdType = c("fdr","pvalue"), maxgap = 500, cutoff = 0.05, minregion = 1) {
    ## set parameters
    thresholdType <- match.arg(thresholdType)
    peakType <- match.arg(peakType)
    if(is.null(threshold)) {
        threshold <- switch(thresholdType,
                            fdr = 0.1,
                            pvalue = 0.05)
    }
    
    splines <- fit@smooths$splines
    index <- fit@smooths$chunkIndex
    if(is.null(smooth)) {
        lables <- colnames(fit@experimentDesign)
        smooth <- paste("s(x)", lables[length(lables)], sep = ":")
    }
    else {
        if(smooth == "") {
            smooth <- "s(x)"
        }
        else {
            smooth <- paste("s(x)", smooth, sep = ":")
        }
    }

    if(is.null(range)) {
        xid <- 1:length(fit@positions)
    }
    else {
        xid <- queryHits(findOverlaps(fit@positions, range))
    }
    
    ## find row index by overlaps and subset position vector
    iid <- queryHits(findOverlaps(index, fit@positions[xid]))

  x <- DataFrame(seqnames = Rle(as.factor(seqnames(fit@positions)[xid])),
                              pos = pos(fit@positions)[xid], id = Rle(iid))
    
    if(peakType == "narrow") {
        futile.logger::flog.info("Calling narrow peaks")
        npeaks <- getExtremes(x, splines, smooth)
        if(nrow(npeaks) == 0) {
          return(data.table())
        }
        peaks <- computePeakSignificance(fit, npeaks, x)
        peaks$id <- NULL
    }
    if(peakType == "broad") {
        futile.logger::flog.info("Calling broad peaks")
        zscore <- computeZscore(fit, xid, smooth)
        bpeaks <- callBroadPeaks(zscore, maxgap, cutoff)
        if(length(bpeaks) == 0) {
          return(data.table())
        }
        peaks <- computeBroadPeakSignificance(fit, bpeaks, smooth)
        peaks <- peaks[width >= minregion,]
    }
    
    if(thresholdType == "pvalue") {        
        signif <- peaks[score >= -log(threshold),]
        signif <- signif[order(score, decreasing = TRUE),]
    }
    if(thresholdType == "fdr") {
        signif <- peaks[fdr <= threshold,]
        signif <- signif[order(fdr),]
    }
    return(signif)
}

#' Get extreme points on a piecewise polynomial function for a single tile
#' This is a function to be applied over in getExtremes. For parameters
#' see getExtremes
#' @noRd
computeTileExtremes <- function(ii, x, splines, smooth, what, useIntercept) {
    intercept <- 0
    
    subx <- data.table::copy(x[id == ii, pos])
    chrom <- data.table::copy(x[id == ii, seqnames])
    sp <- splines[[as.character(ii)]][splines[[as.character(ii)]]$smooth == smooth,]
    if(useIntercept) {
        intercept <- attr(splines[[as.character(ii)]], "intercept")
    }
    
    ## compute basis
    basis <- bspline(subx, sp$knots)
    basis_prime <- bspline(subx, sp$knots, derivative = 1)
    basis_pprime <- bspline(subx, sp$knots, derivative = 2)
    
    ## compute splines
    spline <- computeSpline(basis, na.omit(sp$coefs), intercept = intercept)
    spline_prime <- computeSpline(basis_prime, na.omit(sp$coefs), intercept = intercept)
    spline_pprime <- computeSpline(basis_pprime, na.omit(sp$coefs), intercept = intercept)
    
    ## find roots
    f_prime <- cbind(subx, spline_prime)
    f_pprime <- cbind(subx, spline_pprime)
    roots_prime <- find_root(f_prime)
    roots_pprime <- find_root(f_pprime)
    
    ## classify
    bool_indx <- ((f_pprime[,1] %in% round(roots_prime)) & (f_pprime[,2] < 0))
    if(sum(bool_indx) == 0) {
        extMax <- data.table::data.table(seqnames = factor(), position = integer(), summit = numeric(), type = factor())
    }
    else {
        extMax <- data.table::data.table(seqnames = chrom[bool_indx], position = subx[bool_indx],
                                         summit = exp(spline[bool_indx]), type = "max")
    }
    
    bool_indx <- ((f_pprime[,1] %in% round(roots_prime)) & (f_pprime[,2] > 0))
    if(sum(bool_indx) == 0) {
        extMin <- data.table::data.table(seqnames = factor(), position = integer(), summit = numeric(), type = factor())
    }
    else {
        extMin <- data.table::data.table(seqnames = chrom[bool_indx], position = subx[bool_indx],
                                         summit = exp(spline[bool_indx]), type = "min")
    }
    ext <- rbind(extMax,extMin)
    
    if (what == "all") {
        bool_indx <- ((f_pprime[,1] %in% round(roots_pprime)))
        if(sum(bool_indx) == 0) {
            extInf <- data.table::data.table(seqnames = factor(), position = integer(), summit = numeric(), type = factor())
        }
        else {
            extInf <- data.table::data.table(seqnames = chrom[bool_indx], position = subx[bool_indx],
                                             summit = exp(spline[bool_indx]), type = "inflection")
        }
        ext <- rbind(ext, extInf)
    }
    rm(list = c("subx", "chrom", "sp", "basis", "basis_prime", "basis_pprime", "spline", "spline_prime",
                "spline_pprime", "f_prime", "f_pprime", "roots_prime", "roots_pprime", "bool_indx", "extMin", "extMax"))
    return(ext)
}

#' Get extreme points on a piecewise polynomial function.
#'
#' @param x A numeric vector to be evaluated by the function
#' @param splines The 'smooths' slot from the \code{GenoGAM} object
#' @param smooth The specific smooth name, as specified by the design
#' @param border The border of the peaks defined either by the minimum (default)
#' or the inflection points
#' @return A data.table of all extreme points.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
getExtremes <- function(x, splines, smooth) {q

    ## define variables
    
    intercept <- 0
    useIntercept <- FALSE
    if(smooth == "s(x)") {
        useIntercept <- TRUE
    }
    
    ids <- unique(x$id)

  lambdaFun <- function(params, x, smooth) {
    suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))

    sp <- params[params$smooth == smooth,]
    ii <- attr(params, "id")
    
    subx <- x[x$id == ii, "pos"]
    chrom <- x[x$id == ii, "seqnames"]

    if(useIntercept) {
      intercept <- attr(params, "intercept")
    }
    
    ## compute basis
    basis <- bspline(subx, sp$knots)
    basis_prime <- bspline(subx, sp$knots, derivative = 1)
    basis_pprime <- bspline(subx, sp$knots, derivative = 2)
    
    ## compute splines
    spline <- computeSpline(basis, na.omit(sp$coefs), intercept = intercept)
    spline_prime <- computeSpline(basis_prime, na.omit(sp$coefs))
    spline_pprime <- computeSpline(basis_pprime, na.omit(sp$coefs))
    
    ## find roots
    f_prime <- cbind(subx, spline_prime)
    f_pprime <- cbind(subx, spline_pprime)
    roots_prime <- find_root(f_prime)
    roots_pprime <- find_root(f_pprime)
    
    ## classify
    bool_indx <- ((f_pprime[,1] %in% round(roots_prime)) & (f_pprime[,2] < 0))
    if(sum(bool_indx) == 0) {
      extMax <- data.table::data.table(seqnames = factor(), position = integer(), summit = numeric(), type = factor())
    }
    else {
      extMax <- data.table::data.table(seqnames = as.factor(chrom[bool_indx]), position = as.vector(subx[bool_indx]),
                                       summit = exp(spline[bool_indx]), type = "max")
    }
    
    bool_indx <- ((f_pprime[,1] %in% round(roots_prime)) & (f_pprime[,2] > 0))
    if(sum(bool_indx) == 0) {
      extMin <- data.table::data.table(seqnames = factor(), position = integer(), summit = numeric(), type = factor())
    }
    else {
      extMin <- data.table::data.table(seqnames = as.factor(chrom[bool_indx]), position = as.vector(subx[bool_indx]),
                                       summit = exp(spline[bool_indx]), type = "min")
    }
    ext <- rbind(extMax,extMin)
    
    return(ext)
  }
    
    res <- lapply(splines, lambdaFun, x = x, smooth = smooth)
    
    res <- data.table::rbindlist(res)
    res$type <- as.factor(res$type)
    res <- res[order(position),]
    attr(res, "smooth") <- smooth
    return(res)  
}

#' Compute spline from basis matrix and coefficients
#'
#' @param base A design matrix of size n x k, with n data points evaluated by every base function /eqn{b_k}
#' @param coefs the coefficients of the base functions
#' @param intercept An intercept if necessary (usually not)
#' @return A vector of length n, giving the function values for given data points
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
computeSpline <- function(base, coefs, intercept = 0) {
    base %*% coefs + intercept
}


#' Find roots in a function f.
#'
#' A function to find roots in a function f.
#' 
#' @param f is a data.frame with two columns. The first column is the x or input column.
#' The second is the y or output column
#' @return A numeric vector of x values which have a root
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
find_root <- function(f) {
    
    f_binary <- f[,2] > 0
    dif <- as.logical(abs(diff(f_binary,1)))
    return (f[c(FALSE,dif),1])
}

#' Computing pointwise zscores for given positions
#'
#' @param fit A \code{GenoGAM} object
#' @param index A vector of position indeces in the fit
#' @param smooth The complete name of the smooth, as given
#' in the column names of the 'fits' slot of the \code{GenoGAM} object
#' @return A data.table of positions and their zscores and pvalues
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
computeZscore <- function(fit, index, smooth) {
    
    se_smooth <- paste("se", smooth, sep = ".")
    all <- slot(fit, "fits")[,smooth]
    se_all <- slot(fit, "fits")[,se_smooth]
    
    isInstalled <- require(genefilter, quietly = TRUE)
    if(isInstalled) {
        mu0 <- genefilter::shorth(all, na.rm=TRUE)
    }
    else {
        futile.logger::flog.info("'genefilter' package not installed. Using 'median' to estimate signal average")
        mu0 <- median(all, na.rm=TRUE)
    }
    
    left <- all[all <= mu0]
    right <- abs(left - mu0) + mu0
    new_data <- c(left, right)
    var0 <- mad(new_data, na.rm = TRUE)^2
    zscore <- (all[index] - mu0)/(sqrt(se_all[index]^2 + var0))
        
    res <- data.table::data.table(seqnames = as.factor(seqnames(slot(fit, "positions")[index,])),
                      pos = pos(slot(fit, "positions")[index,]),
                      zscore = zscore)
    res$pval <- -pnorm(-res$zscore, log.p = TRUE)
    return(res)
}


#' Call broad peaks based on pointwise pvalues
#'
#' @param zscore A data.table of positions and their respective zscores and pvalues
#' @param maxgap The maximum gap between two broad peaks, that can be tolerated in order
#' to identify both as part of one broad peak. All broad peaks with distances smaller or
#' equal to the maxgap will be merged.
#' @param cutoff A threshold that is used to identify all significantly high positions,
#' which then make up a broad peak.
#' @return A \code{GPos} object of positions (and regions) that define the broad peaks
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
callBroadPeaks <- function(zscore, maxgap, cutoff) {
    signifPos <- zscore[pval >= -log(cutoff), pos]
    if(length(signifPos) == 0) {
      gp <- GenomicRanges::GPos(GRanges())
    }
    else {
      diffs <- abs(diff(signifPos))
      breaks <- which(diffs > maxgap)
    
      starts <- c(signifPos[1], signifPos[breaks + 1])
      ends <- c(signifPos[breaks], signifPos[length(signifPos)])
      chroms <- c(as.character(zscore[breaks, seqnames]), as.character(zscore[pos == ends[length(ends)], seqnames]))
      gr <- GenomicRanges::GRanges(chroms, IRanges(starts, ends))
      gp <- GenomicRanges::GPos(gr)
      indx <- match(pos(gp), zscore$pos)
      gp$zscore <- zscore[indx, zscore]
      gp$pval <- zscore[indx, pval]
    }
    
    return(gp)
}

#' Computes significance of the narrow peaks
#'
#' @param fit A \code{GenoGAM} object
#' @param peaks A data.table of peaks
#' @param x A data.table with positions and their chunk ids
#' @return A data.table of peaks and their respective information
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
computePeakSignificance <- function(fit, peaks, x) {
  if(length(peaks) == 0) {
    return(peaks)
  }
  smooth <- attr(peaks, "smooth")

  gr <- GRanges(peaks$seqnames, IRanges(peaks$position, peaks$position))
  ov <- findOverlaps(gr, slot(fit, "positions"))
  xid <- subjectHits(ov)
  z <- computeZscore(fit, xid, smooth)
  peaks$zscore <- z$zscore
  
  p <- data.table::copy(peaks[type == "max",])
  v <- data.table::copy(peaks[type == "min",])
  p <- p[order(zscore, decreasing = TRUE),]
  v$zscore <- -v$zscore
  v <- v[order(zscore, decreasing = TRUE),]
  fdr <- sapply(p$zscore, function(y) {
    sum(v$zscore >= y)/sum(p$zscore >= y)
  })
  fdr[is.nan(fdr)] <- 0
  if(nrow(v) > nrow(p)) fdr <- fdr*(nrow(p)/nrow(v)) ## correct for different size of valleys and peaks to avoid FDR > 1
  peaks <- p
  peaks$type <- NULL
  peaks$score <- -pnorm(-peaks$zscore, log.p = TRUE)
  peaks$fdr <- fdr
  attr(peaks, "smooth") <- smooth
  peaks <- xsd(fit, peaks)
  return(peaks)
}

#' Parses data.table of extermums to a data.table of peaks
#'
#' @param peakDF A data.table of exterme points
#' @param x A data.table with positions and their chunk ids
#' @return A data.table of peaks and their respective information
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
parsePeaks <- function(peakDF, x) {
    if(nrow(peakDF) == 0) {
        peaks <- peakDF
        peaks$from <- numeric()
        peaks$to <- numeric()
        peaks$type <- NULL
        return(peaks)
    }
    
    res <- lapply(levels(peakDF$seqnames), function(y) {
        subpeaks <- data.table::copy(peakDF[seqnames == y,])
        orderedPos <- subpeaks[order(position),]
        indx <- which(orderedPos$type == "max")
        peaks <- orderedPos[indx,]
        
        if(orderedPos$type[1] == "max") {
            peaks$start <- c(min(x[seqnames == y, pos], na.rm = TRUE), orderedPos[indx[-1] - 1, position])
        }
        else {
            peaks$start <- orderedPos[indx - 1, position]
        }
        
        if(orderedPos$type[nrow(orderedPos)] == "max") {
            peaks$end <- c(orderedPos[indx[-(length(indx))] + 1, position], min(x[seqnames == y, pos], na.rm = TRUE))
        }
        else {
            peaks$end <- orderedPos[indx + 1, position]
        }
        
        peaks$type <- NULL
        return(peaks)
    })
    res <- data.table::rbindlist(res)
    res <- res[order(zscore, decreasing = TRUE),]
    return(res)
}

#' Computes the variance of the peak locations and a 95\% confidence interval around them
#'
#' @param fit A \code{GenoGAM} object
#' @param peaks An already parsed data.table of peaks
#' @return A data.table of peaks with additional columns of the left and right confidence interval borders
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
xsd <- function(fit, peaks) {
  if(length(peaks) == 0) {
    return(peaks)
  }
  smooth <- attr(peaks, "smooth")
  splines <- slot(fit, "smooths")$splines
  index <- slot(fit, "smooths")$chunkIndex
  v <- slot(fit, "vcov")
  gr <- GRanges(peaks$seqnames, IRanges(start = peaks$position, end = peaks$position))
  id <- subjectHits(findOverlaps(gr, index))
  peaks$id <- id
  peaks <- peaks[order(id),]
  
  res <- sapply(unique(peaks$id), function(ii) {
    p <- peaks[id == ii, position]
    params <- splines[[as.character(ii)]]
    k <- params[params$smooth == smooth, "knots"]
    coefs <- na.omit(params[params$smooth == smooth, "coefs"])
    bsprime <- bspline(p, k, derivative = 1)
    f2 <- bspline(p, k, derivative = 2) %*% coefs
    n <- length(v[[as.character(ii)]])
    
    ## from vector representing the upper triangular vcov matrix
    ## to the full symmetric matrix
    solve_n <- (sqrt(1 + 4*n*2) - 1)/2 ## re-arranged midnight formula to get dimensions of symmetric matrix
    subv <- diag(solve_n)
    subv[upper.tri(subv, diag = TRUE)] <- v[[as.character(ii)]]
    subv <- subv + t(subv) - diag(diag(subv))
    indx <- which(attr(v[[as.character(ii)]], "smooths") == smooth)
    
    betav <- subv[indx, indx]
    xvar <- diag(bsprime %*% betav %*% t(bsprime))/f2^2
    return(xvar)
  })
    
  res <- unlist(res)
  peaks$start <- peaks$position - 1.96 * sqrt(res)
  peaks$end <- peaks$position + 1.96 * sqrt(res)
  peaks$start[which(peaks$start < 1)] <- 1
  for(chrom in seqlevels(index)) {
    indx <- which(peaks[seqnames == chrom, end] > seqlengths(index)[chrom])
    if(length(indx) > 0) {
      pos <- peaks[seqnames == chrom, position][indx]
      peaks[seqnames == chrom & position == pos]$end <- seqlengths(index)[chrom]
    }
  }
  return(peaks)
}

#' Computes the significance of broad peaks
#'
#' @param fit A \code{GenoGAM} object
#' @param peakDF A \code{GPos} object of broad peaks
#' @param smooth The name of the smooth as given in the experiment design
#' @return A data.table of broad peaks with the respective information
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
computeBroadPeakSignificance <- function(fit, peakDF, smooth) {
  if(length(peakDF) == 0) {
    return(data.table())
  }
  bp <- peakDF
  regions <- peakDF@pos_runs
  indx <- queryHits(findOverlaps(slot(fit, "positions"), regions))
  bp$region <- subjectHits(findOverlaps(bp, regions))
  bp$estimate <- exp(slot(fit, "fits")[indx, smooth])
  bp <- data.table::data.table(as.data.frame(bp))
  regions <- data.table::data.table(as.data.frame(regions))
  pv <- bp[, min(p.adjust(exp(-pval), method="hochberg")), by = region]
  fpb <- bp[, mean(estimate), by = region]
  regions$score[pv$region] <- pv$V1
  regions$meanSignal[fpb$region] <- fpb$V1
  regions$fdr = p.adjust(regions$score, method="BH")
  regions$score <- -log(regions$score)
  return(regions)
}
