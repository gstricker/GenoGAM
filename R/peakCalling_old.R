## ======================
## PEAK CALLING
## =====================

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

#' Compute peak significance.
#'
#' A function to compute peak significance.
#'
#' @param fun A function object as returned by the genogam function
#' @param data A list of GenomicTiles on which to evaluate the function
#' @param form A character object, defining a formula or a single spline. See details
#' @param peaks A data.frame of peaks to compute significance for, as returned by .getExtremes
#' @return A data.frame of peaks and additional information
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
.computeScore <- function(fun, data, form, peaks, valleys, intercept = TRUE){
    if(class(data) != "list") data <- list(data)
    if(nrow(peaks) == 0) return(peaks)
    valleys <- valleys[,summit := -summit]

    res <- bplapply(data, function(z) {   
                        require(genoGAM, quietly = TRUE)
                        x <- unique(as.vector(untile(z)$x))
                        chrom <- as.character(unique(seqnames(z)))
                        spline <- fun(x, form, chrom, confidence = TRUE, addIntercept = intercept)
                        spline$sd <- abs(spline$estimate - spline$lower)/1.96
                        spline$position <- x
                        spline$chromosome <- chrom
                        return(data.table(spline))
                    })
    all <- rbindlist(res)
            
    mu0 <- tryCatch({
        shorth(all$estimate, na.rm=TRUE)
    }, error = function(e) {
        median(all$estimate, na.rm=TRUE)
    })

    left <- all[estimate <= mu0,]
    right <- data.table(left)
    right$estimate <- abs(left$estimate - mu0) + mu0
    new_data <- rbind(left,right)
    varmu0 <- mad(new_data$estimate, na.rm = TRUE)^2

    all$zscore <- (all$estimate - mu0)/(sqrt(all$sd^2 + varmu0))
    all$score <- -pnorm(-all$zscore, log = TRUE)
    all$estimate <- NULL

    chromosomes <- unique(all$chromosome)
    indices <- bplapply(chromosomes, function(z) {
        require(genoGAM, quietly = TRUE)
        subdata <- all[chromosome == z,]
        subpeaks <- peaks[chromosome == z,]
        combined <- data.table(merge(subpeaks,subdata,by = "position"))
        combined$chromosome.y <- NULL
        setnames(combined,"chromosome.x","chromosome")
        return(combined)
    })
    allCombined <- rbindlist(indices)
    allCombined[, lower := exp(lower)]
    allCombined[, upper := exp(upper)]

    indices <- bplapply(chromosomes, function(z) {
        require(genoGAM, quietly = TRUE)
        subdata <- all[chromosome == z,]
        subvalleys <- valleys[chromosome == z,]
        combined <- data.table(merge(subvalleys,subdata,by = "position"))
        return(combined)
    })
    combinedValleys <- rbindlist(indices)
    
    allCombined <-  allCombined[order(zscore, decreasing = TRUE),]
    combinedValleys <- combinedValleys[order(zscore, decreasing = TRUE),]
    allCombined$fdr <- sapply(allCombined$zscore, function(y) sum(combinedValleys$zscore > y)/sum(allCombined$zscore > y))*(nrow(allCombined)/nrow(combinedValleys))
    allCombined[is.nan(fdr), fdr:=0]
    
    
    return(allCombined)
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
.callPeaks <- function(fun, data, form, borders = c("min", "inflection"), cutoff = NULL, cutoffType = c("fdr","pvalue"), garbageCollection = FALSE, intercept = TRUE) {
    if(class(data) != "list") data <- list(data)
    cutoffType <- match.arg(cutoffType)
    borders <- match.arg(borders)
    peakBorders <- switch(borders,
                          min = "root",
                          inflection = "scew")

    message("Calling peaks and computing their borders")
    peaks <- bplapply(1:length(data), function(y) {

        require(genoGAM, quietly = TRUE)
        
        ## numWorker <- bpparam()$workers
        ## if(y%%numWorker == length(data)%%numWorker | y == length(data)) {
        ##     progress <- y/length(data)*100
        ##     progressBar <- c("[", rep("*",round(progress)), rep(" ",100 - round(progress)), "]  ", as.character(round(progress,2)) , "%")
        ##     message(paste(progressBar, collapse = ""))
        ##     flush.console()
        ## }
                    
        x <- unique(as.vector(data[[y]]$x))
        chrom <- as.character(unique(seqnames(data[[y]])))
        pos <- .getExtremes(fun, x, form, chrom, what = peakBorders, intercept = intercept)
        peaks <- .parsePeaks(pos, x)
        if(nrow(peaks) != 0) {
            peaks$chromosome <- chrom
            peaks$track <- form
        }
        else {
            peaks$chromosome <- character()
            peaks$track <- character()
        }
        if(garbageCollection) gc()
        return(data.table(peaks))
    })
    
    allPeaks <- rbindlist(peaks)
    uniquePositions <- unique(allPeaks$position)
    allPeaks <- allPeaks[match(uniquePositions, allPeaks$position),]

    message("Calling valleys for FDR computation")
    valleys <- bplapply(1:length(data), function(y) {

        require(genoGAM, quietly = TRUE)
        
        ## numWorker <- bpparam()$workers
        ## if(y%%numWorker == length(data)%%numWorker | y == length(data)) {
        ##     progress <- y/length(data)*100
        ##     progressBar <- c("[", rep("*",round(progress)), rep(" ",100 - round(progress)), "]  ", as.character(round(progress,2)) , "%")
        ##     message(paste(progressBar, collapse = ""))
        ##     flush.console()
        ## }
                    
        x <- unique(as.vector(data[[y]]$x))
        chrom <- as.character(unique(seqnames(data[[y]])))
        pos <- .getExtremes(fun, x, form, chrom, what = "min", intercept = intercept)
        pos$chromosome <- chrom
        
        if(garbageCollection) gc()  
        return(data.table(pos))
    })

    allValleys <- rbindlist(valleys)
    uniquePositions <- unique(allValleys$position)
    allValleys <- allValleys[match(uniquePositions, allValleys$position),]

    message("Computing significance scores")
    signif <- .computeScore(fun, data, form, allPeaks, allValleys, intercept = intercept)
    if(nrow(signif) == 0) return(signif)
    if(!is.null(cutoff)) {
        if(cutoffType == "pvalue") signif <- signif[score >= -log10(cutoff),]
        if(cutoffType == "fdr") signif <- signif[fdr <= cutoff,]
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
.getExtremes <- function(fun, x, form, chrom, what = c("max", "min", "inflection", "all", "root", "scew"), intercept = TRUE) {
    
    splines_first_deriv <- 0
    splines_second_deriv <- 0
    spline <- fun(x, form, chrom, addIntercept = intercept)
    splines_first_deriv <- fun(x, form, chrom, deriv = 1)
    splines_second_deriv <- fun(x, form, chrom, deriv = 2)
    
    f <- cbind(x, splines_first_deriv$estimate)
    f_prime<- cbind(x, splines_second_deriv$estimate)
    
    roots <- .find_root(f)

    what <- match.arg(what)
    
    if(what == "inflection" | what == "all" | what == "scew") inflection <- .find_root(f_prime)

    if (what == "max") {
        bool_indx <- ((f_prime[,1] %in% round(roots)) & (f_prime[,2] < 0))
        if(sum(bool_indx) == 0) ext <- data.table(position = numeric(), summit = numeric(), type = character())
        else ext <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "max")
    }
    if (what == "min") {
        bool_indx <- ((f_prime[,1] %in% round(roots)) & (f_prime[,2] > 0))
        if(sum(bool_indx) == 0) ext <- data.table(position = numeric(), summit = numeric(), type = character())
        else ext <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "min")
    }
    if (what == "inflection") {
        bool_indx <- f_prime[,1] %in% round(inflection)
        if(sum(bool_indx) == 0) ext <- data.table(position = numeric(), summit = numeric(), type = character())
        else ext <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "inflection")
    }
    if (what == "all") {
        bool_indx <- ((f_prime[,1] %in% round(roots)) & (f_prime[,2] < 0))
        if(sum(bool_indx) == 0) extMax <- data.table(position = numeric(), summit = numeric(), type = character())
        else extMax <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "max")
        bool_indx <- ((f_prime[,1] %in% round(roots)) & (f_prime[,2] > 0))
        if(sum(bool_indx) == 0) extMin <- data.table(position = numeric(), summit = numeric(), type = character())
        else extMin <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "min")
        bool_indx <- f_prime[,1] %in% round(inflection)
        if(sum(bool_indx) == 0) extInf <- data.table(position = numeric(), summit = numeric(), type = character())
        else extInf <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "inflection")
        ext <- rbind(extMax,extMin,extInf)
    }
    if (what == "root") {
        bool_indx <- ((f_prime[,1] %in% round(roots)) & (f_prime[,2] < 0))
        if(sum(bool_indx) == 0) extMax <- data.table(position = numeric(), summit = numeric(), type = character())
        else extMax <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "max")
        bool_indx <- ((f_prime[,1] %in% round(roots)) & (f_prime[,2] > 0))
        if(sum(bool_indx) == 0) extMin <- data.table(position = numeric(), summit = numeric(), type = character())
        else extMin <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "min")
        ext <- rbind(extMax,extMin)
    }
    if (what == "scew") {
        bool_indx <- ((f_prime[,1] %in% round(roots)) & (f_prime[,2] < 0))
        if(sum(bool_indx) == 0) extMax <- data.table(position = numeric(), summit = numeric(), type = character())
        else extMax <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "max")
        bool_indx <- f_prime[,1] %in% round(inflection)
        if(sum(bool_indx) == 0) extInf <- data.table(position = numeric(), summit = numeric(), type = character())
        else extInf <- data.table(position = x[bool_indx], summit = exp(spline$estimate[bool_indx]), type = "inflection")
        ext <- rbind(extMax,extInf)
    }
        
    ext <- ext[order(ext$pos),]
    return(ext)  
}

#' Add peak bordes.
#'
#' Reshape data.table of extreme points to a peak data.table containing summit and borders
#' @param peakDF A data.table of extreme points as created by .getExtremes
#' @param x A vector of integer positions on which the extremes were called
#' @return A data.table of peaks including borders
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
.parsePeaks <- function(peakDF, x) {
    if(nrow(peakDF) == 0) {
        peaks <- peakDF
        peaks$from <- numeric()
        peaks$to <- numeric()
        peaks$type <- NULL
        return(peaks)
    }

    orderedPos <- peakDF[order(position),]
    indx <- which(orderedPos$type == "max")
    peaks <- orderedPos[indx,]

    if(orderedPos$type[1] == "max") peaks$from <- c(min(x, na.rm = TRUE), orderedPos[indx[-1] - 1, position])
    else peaks$from <- orderedPos[indx - 1, position]

    if(orderedPos$type[nrow(orderedPos)] == "max") peaks$to <- c(orderedPos[indx[-(length(indx))] + 1, position], max(x, na.rm = TRUE))
    else peaks$to <- orderedPos[indx + 1, position]

    peaks$type <- NULL
    return(peaks)
}
