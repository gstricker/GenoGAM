#' The pot function for a GenoGAM object
#'
#' This functions plots the fit of a given region and optionally the read counts from the GenoGAMDataSet object
#'
#' @param x A \code{GenoGAM} object
#' @param ggd A \code{GenoGAMDataSet} object to plot raw counts
#' @param ranges A \code{GRanges} object specifying a particular region
#' @param seqnames A chromosome name. Together with start and end it is an
#' alternative way ob selecting a region
#' @param start The start of a region
#' @param end The end of a region
#' @param scale Logical, should all tracks be scaled to the same y-axis?
#' @param cap If FALSE deactivates the cap that prevents to accidentaly plot a too larger
#' area. The default cap is 1Mbp.
#' @param log Should log values be plotted on the y-axis 
#' @param ... Additional parameters that will be passed to the basic plot routine
#' @return A plot of all tracks either using the ggplot2 or the base R framework
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export

plotGenoGAM <- function(x, ggd = NULL, ranges = NULL, seqnames = NULL,
                         start = NULL, end = NULL, scale = TRUE, cap = TRUE,
                         log = FALSE,...) {
    
    ## determine what type of object we are dealing with
    is_hdf5 <- is.HDF5(x)
    if(is(x, "GenoGAMList")) {
        is_split <- TRUE
    }
    else {
        if(is(x, "GenoGAM")) {
            is_split <- FALSE
        }
        
        else {
            stop("Wrong class submitted")
        }
    }
    
    ## Set cap if not FALSE
    if(!cap) {
        cap <- Inf
    }
    else {
        ## 1Mio points should be fine
        cap <- 1e6
    }
    
    if(is.null(seqnames) & is.null(start) & is.null(end) & is.null(ranges)) {
        stop("Either 'seqnames', and/or 'start' and/or 'end' or 'ranges' has to be provided")
    }
    else {
        ## take care of split rowRanges and different ways
        ## of supplying the arguments
        loc <- rowRanges(x)
        
        ## ranges is preferred to seqnames, start, end
        if(!is.null(ranges)) {
            if(length(ranges) > 1) {
                warning("Only one range can be plotted. Taking the first one.")
            }
            if(is_split) {
                seqnames <- as.character(GenomicRanges::seqnames(ranges)[1])
                loc <- loc[[seqnames]]
            }
            ov <- IRanges::findOverlaps(loc, ranges[1])
            indx <- S4Vectors::queryHits(ov)
        }
        ## if ranges are not provided revert to seqnames and optional start and end
        else {
            if(is.null(seqnames)) {
                stop("Neither 'ranges' nor 'seqnames' are provided. At least one is needed.")
            }
            if(is.null(start)) {
                warning("'start' is not provided. Setting it to 1.")
                start <- 1
            }
            if(is.null(end)) {
                warning(paste("'end' is not provided. Setting it to the last position of chromosome", seqnames))
                end <- Inf
            }
            if(is_split) {
                loc <- loc[[seqnames]]
            }
            indx <- which(GenomicRanges::seqnames(loc) == seqnames &
                          GenomicRanges::pos(loc) >= start & GenomicRanges::pos(loc) <= end)
        }
    }
    if(length(indx) > cap) {
        stop("The fit is too big. Plotting is not advised, please provide a smaller region.")
    }
  
    ## assemble the fit
    pos <- GenomicRanges::pos(loc)[indx]
    if(is_split){
        y <- fits(x)[[seqnames]][indx,]
    }
    else {
        y <- fits(x)[indx,]
    }

    ## assemble the standard errors
    if(is_split){
        se <- se(x)[[seqnames]][indx,]
    }
    else {
        se <- se(x)[indx,]
    }

    ## take care of raw data
    inputData <- NULL ## the raw data, not present by default
    if(!is.null(ggd)) {
        if(is_split){
            inputData <- assay(ggd)[[seqnames]][indx,]
        }
        else {
            inputData <- assay(ggd)[indx,]
        }
        ## normalize by size factors
        if(is_hdf5) {
            aux <- lapply(1:ncol(inputData), function(y) {
                inputData[,1]/exp(sizeFactors(ggd)[y])
            })
            names(aux) <- colnames(inputData)
            inputData <- DataFrame(aux)
        }
        else {
            inputData[] <- lapply(colnames(inputData), function(y) {
                inputData[[y]]/exp(sizeFactors(ggd)[y])
            })
        }
    }
    title <- as.character(.extractGR(loc[indx]))

    ## if(require(Gviz)) {
    ##   plot_gviz()
    ## }
    ## if(requireNamespace("ggplot2", quietly = TRUE) & requireNamespace("grid", quietly = TRUE)) {
    ## plot_ggplot2(x = pos, y = y, se = se, inputData = inputData, scale = scale, title = title)
    ## }
    ## else {
    plot_base(x = pos, y = y, se = se, inputData = inputData,
              scale = scale, title = title, log, ...)
    ## }
}

#' @noRd
plot_base <- function(x, y, se, inputData = NULL, scale = TRUE,
                      title = "", log = FALSE, ...) {

    numTracks <- ncol(y)
    if(!is.null(inputData)) {
        numRaw <- ncol(inputData)
    }

    ## compute confidence interval
    intervals <- vector("list", numTracks)
    names(intervals) <- colnames(y)
    for(name in colnames(y)) {
        if(log) {
            intervals[[name]]$upper <- y[,name] + 2*se[,name]
            intervals[[name]]$lower <- y[,name] - 2*se[,name]
        }
        else {
            intervals[[name]]$upper <- exp(y[,name] + 2*se[,name])
            intervals[[name]]$lower <- exp(y[,name] - 2*se[,name])
        }
    }

    ## compute ylims  
    ymin <- sapply(intervals, function(x) min(x$lower))
    ymax <- sapply(intervals, function(x) max(x$upper))
    if(!is.null(inputData)) {
        input_min <- sapply(inputData, min)
        input_max <- sapply(inputData, max)
    }
    
    if(scale) {
        fit_ylim <- c(min(ymin), max(ymax))
        if(!is.null(inputData)) {
            input_ylim <- c(min(input_min), max(input_max))   
        }
    }
    else {
        fit_ylim <- list(ymin, ymax)
        names(fit_ylim) <- c("min", "max")
        if(!is.null(inputData)) {
            input_ylim <- list(input_min, input_max)
            names(input_ylim) <- c("min", "max")
        }
    }

    
    ## set penalty against scientific numbers
    opt <- options("scipen" = 100000)
    
    ## plot raw counts
    if(!is.null(inputData)) {
        par(mfrow = c(numRaw, 1), oma = c(0, 0, 2, 0))
        for (ii in 1:numRaw) {
            ## write xlab only at the bottom plot
            xlab = ""
            if(ii == numRaw) {
                xlab = "Genomic Position"
            }
            if(!scale){
                ylim = c(input_ylim[['min']][ii], input_ylim[['max']][ii])
            }
            else {
                ylim <- input_ylim
            }
            plot(x, inputData[,ii], type = "p", col = "#73737330", pch = 19,
                 ylim = ylim, xlab = xlab, ylab = colnames(inputData)[ii])
        }
        title(main = title,outer=TRUE)
    }

    ## plot fits and confidence interval
    if(!is.null(inputData)) {
        x11()
    }
    par(mfrow = c(numTracks, 1), oma = c(0, 0, 2, 0))
    for (ii in 1:numTracks) {
        ## write xlab only at the bottom plot
        xlab = ""
        if(ii == numTracks) {
            xlab = "Genomic Position"
        }
        ## scale accordingly
        if(!scale){
            ylim = c(fit_ylim[['min']][ii], fit_ylim[['max']][ii])
        }
        else {
            ylim <- fit_ylim
        }
        ## plot log or normal fit
        if(log) {
            plot(x, y[,ii], type = "l", col = "black", ylim = ylim,
                 xlab = xlab, ylab = colnames(y)[ii], ...)
        }
        else {
            plot(x, exp(y[,ii]), type = "l", col = "black", ylim = ylim,
                 xlab = xlab, ylab = colnames(y)[ii], ...)
        }
        ## plot confidence interval
        lines(x, intervals[[ii]]$lower, lty = "dotted", col = 'grey')
        lines(x, intervals[[ii]]$upper, lty = "dotted", col = 'grey')
        ## plot 0 (log) or 1 (normal) horizontal line
        if(log) {
            abline(h = 0, col = "red")
        }
        else {
            abline(h = 1, col = "red")
        }
    }
    title(main = title,outer=TRUE)
    ## unset penalty for scientific numbers
    options(opt)
}

## #' @noRd
## plot_ggplot2 <- function(x, y, inputData = NULL, scale = TRUE, title = "") {
##   numTracks <- ncol(y)/2

##   seCols <- grep("se", names(y))
##   pCols <- grep("pvalue", names(y))
##   sCols <- (1:ncol(y))[-c(seCols, pCols)]
  
##   ## compute confidence interval
##   for(ii in seCols) {
##     secolname <- names(y)[ii]
##     scolname <- strsplit(secolname, "\\.")[[c(1,2)]]
    
##     uname <- paste("upper", scolname, sep = ".")
##     lname <- paste("lower", scolname, sep = ".")
##     y[[uname]] <- y[[scolname]] + 1.96*y[[secolname]]
##     y[[lname]] <- y[[scolname]] - 1.96*y[[secolname]]
##   }
  
##   ## compute ylim
##   fit_ylim <- NULL
##   input_ylim <- NULL

##   if(scale) {
##     upperCols <- grep("upper.s\\(x\\)", names(y))
##     lowerCols <- grep("lower.s\\(x\\)", names(y))
##     fit_ylim <- c(min(y[,lowerCols]), max(y[,upperCols]))
##     if(!is.null(inputData)) {
##       input_ylim <- c(min(sapply(inputData, min)), max(max(sapply(inputData, max))))
##     }
##   }
    
##   ## put all in one data.frame for ggplot
##   if(!is.null(inputData)) {
##     numTracks <- numTracks + ncol(inputData)
##     y <- cbind(y, as.data.frame(inputData), x)
##   }
##   else {
##     y <- cbind(y, x)
##   }

##   ## create List where to store plots
##   plotList <- vector("list", numTracks)
##   idx <- 1

##   ## plot raw data to List
##   if(!is.null(inputData)) {
##     for (ii in 1:ncol(inputData)) {
##       inputName <- names(inputData)[ii]
##       if(is.null(input_ylim)) {
##         input_ylim <- range(y[[inputName]])
##       }
##       assign(paste0("yinput", idx), y[[inputName]])
##       plotList[[idx]] <- ggplot2::ggplot(y, ggplot2::aes(x = x, y = get(paste0("yinput", idx)))) + 
##         ggplot2::geom_point(color = "#73737330") + ggplot2::ylim(input_ylim) + 
##       ggplot2::ylab(inputName) + ggplot2::xlab("")
##       idx <- idx + 1
##     }
##   }

##   ## plot fits to List
##   for(ii in sCols) {
##     track <- names(y)[ii]

##     ## print 'genomic position' only on the last plot
##     if(ii == sCols[length(sCols)]) {
##       xlab <- "genomic position"
##     }
##     else {
##       xlab <- ""
##     }
##     ## if no ylim provided make sure that confidence intervalls are included
##     if(!scale) {
##       fit_ylim <- c(min(y[[paste("lower", track, sep = ".")]]),
##                 max(y[[paste("upper", track, sep = ".")]]))
##     }
    
##     assign(paste0("yinput", idx), y[[track]])
##     assign(paste0("lower", idx), y[[paste("lower", track, sep = ".")]]) 
##     assign(paste0("upper", idx), y[[paste("upper", track, sep = ".")]]) 
##     plotList[[idx]] <- ggplot2::ggplot(y, ggplot2::aes(x, get(paste0("yinput", idx)))) + 
##       ggplot2::geom_ribbon(ggplot2::aes(ymin = get(paste0("lower", idx)),
##                                         ymax = get(paste0("upper", idx))), fill = "grey70") + 
##       ggplot2::geom_line() + ggplot2::ylim(fit_ylim) + ggplot2::ylab(track) + 
##       ggplot2::xlab(xlab) + ggplot2::geom_hline(yintercept = 0, colour = "red")
##     idx <- idx + 1
##   }

##   ## arrange plots in grid
##   grid::grid.newpage()
##   grid::pushViewport(grid::viewport(layout = grid::grid.layout(numTracks + 1, 1, heights = grid::unit(c(0.5, rep(10/numTracks, numTracks)), "null"))))
##   grid::grid.text(title, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
##   for(idx in 1:length(plotList)) {
##     print(plotList[[idx]], vp = grid::viewport(layout.pos.row = idx + 1, layout.pos.col = 1))         
##   }
## }
