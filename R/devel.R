## ##################
## ## mechanism to store knots and coefficients (see smooths slot in GenoGAM-class)
## ###################

## ## load libraries if necessary
library(GenoGAM)

## ## test data
dir <- "/s/project/coreProm/data/sacCer2"
bppk <- 100
chunksize <- bppk*50
ov <- bppk*10

## read data
ggd <- GenoGAMDataSet("config.txt", chunkSize = chunksize, overhangSize = ov,
                      design = ~ s(x) + s(x, by = tfiib), directory = dir)

gr <- GRanges("chrI", IRanges(100001, 150000))
subggd <- ggd[gr]
fit <- genogam(subggd, bpknots = bppk, lambda = 10, family = mgcv::nb(theta = 2))

range <- GRanges(c("chrI", "chrI"), IRanges(c(105001, 135001), c(122000, 141000)))

#### peak position confidence
n <- 2000
y <- sin(seq(-6,6, length.out = n)) + 10
mu <- y + rnbinom(n, size = 10, mu = y/5)
x <- 1:n
mu <- y + runif(n, min = 1, max = 10)

mod <- mgcv::gam(mu ~ s(x, bs = "ps", sp = 0))
v <- mgcv::vcov.gam(mod)
splines <- extractSplines(mod)
track <- as.vector(splines::spline.des(splines$knots, x, 4, rep(0,length(x)))$design %*% na.omit(splines$coefs))
f1 <- as.vector(splines::spline.des(splines$knots, x, 4, rep(1,length(x)))$design %*% na.omit(splines$coefs))
f <- cbind(x, f1)
roots <- .find_root(f)
bsprime <- splines::spline.des(splines$knots, roots[1], 4, 1)$design
f2 <- as.vector(splines::spline.des(splines$knots, roots[1], 4, 2)$design %*% na.omit(splines$coefs))
xvar <- as.numeric(bsprime %*% v %*% t(bsprime))/f2^2
conf <- c(roots[1] - 1.96 * sqrt(xvar), roots[1] + 1.96 * sqrt(xvar))
plot(x, track, type = "l")
abline(v = roots[c(1,3)], col = "red")
abline(v = conf, col = "blue")

#############################################
#############################################
## peak calling

## compute spline from basis matrix and coefficients
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
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
.find_root <- function(f) {

    f_binary <- f[,2] > 0
    dif <- as.logical(abs(diff(f_binary,1)))
    return (f[c(FALSE,dif),1])
}

#' Call peaks on a piecewise polynomial function.
#'
#' This function looks for peaks on a piecewise polynomial function given the data and
#' a formula of the splines.
#'
#' @param fun A function object as returned by the genogam function
#' @param data A list of GenomicTiles on which to evaluate the function
#' @param form A character object, defining a formula or a single spline. See details
#' @param borders Should peak borders be taken as the minima or inflection points around peaks.
#' @param cutoff The p-value at which to cut off
#' @return A data.frame of peaks and additional data
#' @details The function is preferably used within the genogam function and is made available
#' for the purpose of modularity of the analysis. The input arguments are very specificly relying
#' on the prioir function. Like fun on genogam and data on makeGenomicTiles. The 'form' argument
#' is a character string that serves the different experimental designs. E.g Calling peaks on a single track,
#' let's say, the Input, would be just form = "input". Calling peaks in a more complicated design might be
#' form = "PolII + DifferenceTrack" to call peaks on the modified PollII track.#' 
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
callPeaks <- function(fit, smooth = NULL, range = NULL, peakType = c("narrow", "broad"), border = c("min", "inflection"), threshold = NULL, thresholdType = c("fdr","pvalue"), maxgap = 500, cutoff = 0.05) {
    ## set parameters
    thresholdType <- match.arg(thresholdType)
    border <- match.arg(border)
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
        smooth <- lables[length(lables)]
    }

    if(is.null(range)) {
        xid <- 1:length(fit@positions)
    }
    else {
        xid <- queryHits(findOverlaps(fit@positions, range))
    }
    
    ## find row index by overlaps and subset position vector
    iid <- queryHits(findOverlaps(index, fit@positions[xid]))
    x <- data.table::data.table(seqnames = seqnames(fit@positions)[xid],
                    pos = pos(fit@positions)[xid], id = iid)

    if(peakType == "narrow") {
        futile.logger::flog.info("Calling narrow peaks")
        npeaks <- getExtremes(x, splines, smooth, border)
        peaks <- computePeakSignificance(fit, npeaks, x)
    }
    if(peakType == "broad") {
        futile.logger::flog.info("Calling broad peaks")
        zscore <- computeZscore(fit, xid, smooth)
        bpeaks <- callBroadPeak(zscore, maxgap, cutoff)
        peaks <- computeBroadPeakSignificance(bpeaks)
    }

    if(thresholdType == "pvalue") {        
        signif <- peaks[pvalue >= -log(threshold),]
    }
    if(cutoffType == "fdr") {
        signif <- peaks[fdr <= threshold,]
    }
    return(signif)
}

#' Get extreme points on a piecewise polynomial function.
#'
#' A function to get extreme points on a piecewise polynomial function
#'
#' @param fun A function object as returned by the genogam function
#' @param x A numeric vector to evaluate by the function
#' @param form A character object, defining a formula or a single spline. See details
#' @param chrom The chromosome name where the input of the 'x' argument came from
#' @param what The type of extreme points. Either 'max', 'min', 'inflection',
#' 'root' ('min' + 'max'), 'scew' ('max' + 'inflection')  or 'all'.
#' @return A data.frame of the position and height of the extreme points.
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
#' @noRd
getExtremes <- function(x, splines, smooth, border = c("min", "inflection")) {

    ## define variables
    what <- switch(border,
                   min = "root",
                   inflection = "all")

    if(smooth == "") {
        colname <- "s(x)"
    }
    else {
        colname <- paste("s(x)", smooth, sep = ":")
    }
    
    intercept <- 0
    ids <- unique(x$id)
    
    res <- BiocParallel::bplapply(ids, function(ii) {
        require(GenoGAM, quietly = TRUE)
        
        subx <- x[x$id == ii, pos]
        chrom <- x[x$id == ii, seqnames]
        sp <- splines[[as.character(ii)]][splines[[as.character(ii)]]$smooth == colname,]
        if(colname == "s(x)") {
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
        roots_prime <- .find_root(f_prime)
        roots_pprime <- .find_root(f_pprime)

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
        return(ext)
    }) 

    res <- data.table::rbindlist(res)
    res$type <- as.factor(res$type)
    res <- res[order(res$position),]
    attr(res, "smooth") <- smooth
    return(res)  
}

## computing pointwise zscore for given rowindex
computeZscore <- function(fit, index, smooth) {

    se_smooth <- paste("se.s(x)", smooth, sep = ":")
    smooth <- paste("s(x)", smooth, sep = ":")
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

    res <- data.frame(seqnames = seqnames(slot(fit, "positions")[index,]),
                      pos = pos(slot(fit, "positions")[index,]),
                      zscore = (all[index] - mu0)/(sqrt(se_all[index]^2 + var0)))
    res$pval <- -pnorm(-res$zscore, log.p = TRUE)
    return(res)
}

## call broad peaks based on pointwise pvalue
callBroadPeak <- function(zscore, maxgap, cutoff) {
    signifPos <- zscore[zscore$pval >= -log(cutoff), "pos"]
    diffs <- abs(diff(signifPos))
    breaks <- which(diffs > maxgap)

    starts <- c(signifPos[1], signifPos[breaks + 1])
    ends <- c(signifPos[breaks], signifPos[length(signifPos)])
    chroms <- c(as.character(zscore[breaks, "seqnames"]), as.character(zscore[zscore$pos == ends[length(ends)], "seqnames"]))
    gr <- GenomicRanges::GRanges(chroms, IRanges(starts, ends))
    gp <- GenomicRanges::GPos(gr)
    indx <- match(pos(gp), zscore$pos)
    gp$zscore <- zscore[indx, "zscore"]
    gp$pval <- zscore[indx, "pval"]
    
    return(gp)
}

computePeakSignificance <- function(fit, peaks, x) {
    xid <- match(peaks$position, pos(slot(fit, "positions")))
    z <- computeZscore(fit, xid, attr(peaks, "smooth"))
    peaks$zscore <- z$zscore

    p <- peaks[type == "max",]
    v <- peaks[type == "min",]
    v$zscore <- -v$zscore
    p <- p[order(zscore, decreasing = TRUE),]
    v <- v[order(zscore, decreasing = TRUE),]
    fdr <- sapply(p$zscore, function(y) {
        sum(v$zscore >= y)/(sum(v$zscore >= y) + sum(p$zscore >= y))
    })
    fdr[is.nan(fdr)] <- 0
    peaks <- parsePeaks(peaks, x)
    peaks$pvalue <- -pnorm(-peaks$zscore, log.p = TRUE)
    peaks$fdr <- fdr
    return(peaks)
}

parsePeaks <- function(peakDF, x) {
    if(nrow(peakDF) == 0) {
        peaks <- peakDF
        peaks$from <- numeric()
        peaks$to <- numeric()
        peaks$type <- NULL
        return(peaks)
    }

    if("inflection" %in% levels(peakDF$type)) {
        peakDF <- peakDF[type != "min",]
    }

    res <- lapply(levels(peakDF$seqnames), function(y) {
        subpeaks <- peakDF[seqnames == y,]
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

computeBroadPeakSignificance <- function(peakDF) {
    bp <- peakDF
    regions <- peakDF@pos_runs
    bp$region <- subjectHits(findOverlaps(bp, regions))
    bp <- data.table::data.table(as.data.frame(bp))
    regions <- data.table::data.table(as.data.frame(regions))
    pv <- bp[, min(p.adjust(exp(-pval), method="hochberg")), by = region]
    regions$pvalue[pv$region] <- pv$V1
    regions$fdr = p.adjust(regions$pvalue, method="BH")
    regions$pvalue <- -log(regions$pvalue)
    return(regions)
}

partial_xsd <- function(mod, x) {
    v <- mgcv::vcov.gam(mod)
    num_var <- length(mod$smooth)
    res <- data.table::data.table(matrix(0, length(x), num_var))
    for(jj in 1:num_var) {
        knots <- mod$smooth[[jj]]$knots
        first <- mod$smooth[[jj]]$first.para
        last <- mod$smooth[[jj]]$last.para
        label <- mod$smooth[[jj]]$label
        label <- paste("pxsd", gsub("pos", "x", label), sep = ".")
        if(jj == 1) {
            first <- 1
        }
        ans <- sapply(x, function(y) {
            bsprime <- bspline(y, knots, derivative = 1)
            val <- as.numeric(bsprime %*% v[first:last, first:last] %*% t(bsprime))
            return(val)
        })
        setnames(res, names(res)[jj], label)
        res[[label]] <- ans
    }
    return(res)
}
