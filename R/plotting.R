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
#' @param ... Additional parameters that will be passed to the basic plot routine
#' @return A plot of all tracks either using the ggplot2 or the base R framework
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export

plot.GenoGAM <- function(x, ggd = NULL, ranges = NULL, seqnames = NULL,
                      start = NULL, end = NULL, scale = TRUE, ...) {
  ## Cap for too long regions
  cap <- 1e5
  loc <- rowRanges(x)
  indx <- 1:length(loc)

  if(is.null(seqnames) & is.null(start) & is.null(end) & is.null(ranges)) {
    sub <- x
    if(length(rowRanges(fit)) > cap) {
      stop("The entire fit is too big. Plotting is not advised, please provide a smaller region.")
    }
  }
  else {
    
    if(!is.null(ranges)) {
      ov <- findOverlaps(loc, ranges)
      indx <- queryHits(ov)
    }
    else {
      if(is.null(start)) {
        start <- 1
      }
      if(is.null(end)) {
        end <- Inf
      }
      if(is.null(seqnames)){
        indx <- which(pos(loc) >= start & pos(loc) <= end)
      }
      else {
        indx <- which(seqnames(loc) == seqnames & pos(loc) >= start & pos(loc) <= end)
      }
    }
  }
  if(length(indx) > cap) {
    stop("The fit is too big. Plotting is not advised, please provide a smaller region.")
  }
  
  ## assemble
  pos <- pos(loc)[indx]
  y <- getFits(x)[indx,]
  inputData <- NULL ## the raw data, not present by default
  if(!is.null(ggd)) {
    inputData <- assay(ggd)[indx,]
    inputData[] <- lapply(names(inputData), function(y) inputData[[y]]/exp(sizeFactors(ggd)[y]))
  }
  title <- as.character(loc[indx]@pos_runs)

  ## if(require(Gviz)) {
  ##   plot_gviz()
  ## }
  if(requireNamespace("ggplot2", quietly = TRUE) & requireNamespace("grid", quietly = TRUE)) {
    plot_ggplot2(x = pos, y = y, inputData = inputData, scale = scale, title = title)
  }
  else {
    plot_base(x = pos, y = y, inputData = inputData, scale = scale, title = title, ...)
  }
}

#' @noRd
plot_ggplot2 <- function(x, y, inputData = NULL, scale = TRUE, title = "") {
  numTracks <- ncol(y)/2

  seCols <- grep("se", names(y))
  pCols <- grep("pvalue", names(y))
  sCols <- (1:ncol(y))[-c(seCols, pCols)]
  
  ## compute confidence interval
  for(ii in seCols) {
    secolname <- names(y)[ii]
    scolname <- strsplit(secolname, "\\.")[[c(1,2)]]
    
    uname <- paste("upper", scolname, sep = ".")
    lname <- paste("lower", scolname, sep = ".")
    y[[uname]] <- y[[scolname]] + 1.96*y[[secolname]]
    y[[lname]] <- y[[scolname]] - 1.96*y[[secolname]]
  }
  
  ## compute ylim
  fit_ylim <- NULL
  input_ylim <- NULL

  if(scale) {
    upperCols <- grep("upper.s\\(x\\)", names(y))
    lowerCols <- grep("lower.s\\(x\\)", names(y))
    fit_ylim <- c(min(y[,lowerCols]), max(y[,upperCols]))
    if(!is.null(inputData)) {
      input_ylim <- c(min(sapply(inputData, min)), max(max(sapply(inputData, max))))
    }
  }
    
  ## put all in one data.frame for ggplot
  if(!is.null(inputData)) {
    numTracks <- numTracks + ncol(inputData)
    y <- cbind(y, as.data.frame(inputData), x)
  }
  else {
    y <- cbind(y, x)
  }

  ## create List where to store plots
  plotList <- vector("list", numTracks)
  idx <- 1

  ## plot raw data to List
  if(!is.null(inputData)) {
    for (ii in 1:ncol(inputData)) {
      inputName <- names(inputData)[ii]
      if(is.null(input_ylim)) {
        input_ylim <- range(y[[inputName]])
      }
      assign(paste0("yinput", idx), y[[inputName]])
      plotList[[idx]] <- ggplot2::ggplot(y, ggplot2::aes(x = x, y = get(paste0("yinput", idx)))) + 
        ggplot2::geom_point(color = "#73737330") + ggplot2::ylim(input_ylim) + 
      ggplot2::ylab(inputName) + ggplot2::xlab("")
      idx <- idx + 1
    }
  }

  ## plot fits to List
  for(ii in sCols) {
    track <- names(y)[ii]

    ## print 'genomic position' only on the last plot
    if(ii == sCols[length(sCols)]) {
      xlab <- "genomic position"
    }
    else {
      xlab <- ""
    }
    ## if no ylim provided make sure that confidence intervalls are included
    if(!scale) {
      fit_ylim <- c(min(y[[paste("lower", track, sep = ".")]]),
                max(y[[paste("upper", track, sep = ".")]]))
    }
    
    assign(paste0("yinput", idx), y[[track]])
    assign(paste0("lower", idx), y[[paste("lower", track, sep = ".")]]) 
    assign(paste0("upper", idx), y[[paste("upper", track, sep = ".")]]) 
    plotList[[idx]] <- ggplot2::ggplot(y, ggplot2::aes(x, get(paste0("yinput", idx)))) + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = get(paste0("lower", idx)),
                                        ymax = get(paste0("upper", idx))), fill = "grey70") + 
      ggplot2::geom_line() + ggplot2::ylim(fit_ylim) + ggplot2::ylab(track) + 
      ggplot2::xlab(xlab) + ggplot2::geom_hline(yintercept = 0, colour = "red")
    idx <- idx + 1
  }

  ## arrange plots in grid
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(numTracks + 1, 1, heights = grid::unit(c(0.5, rep(10/numTracks, numTracks)), "null"))))
  grid::grid.text(title, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  for(idx in 1:length(plotList)) {
    print(plotList[[idx]], vp = grid::viewport(layout.pos.row = idx + 1, layout.pos.col = 1))         
  }
}

#' @noRd
plot_base <- function(x, y, inputData = NULL, scale = TRUE, title = "", ...) {

  numTracks <- ncol(y)/2

  seCols <- grep("se", names(y))
  sCols <- (1:ncol(y))[-seCols]
  
  ## compute confidence interval
  for(ii in seCols) {
    secolname <- names(y)[ii]
    scolname <- strsplit(secolname, "\\.")[[c(1,2)]]
    
    uname <- paste("upper", scolname, sep = ".")
    lname <- paste("lower", scolname, sep = ".")
    y[[uname]] <- y[[scolname]] + 1.96*y[[secolname]]
    y[[lname]] <- y[[scolname]] - 1.96*y[[secolname]]
  }

  ## compute ylims
  fit_ylim <- NULL
  input_ylim <- NULL
  
  if(scale) {
    upperCols <- grep("upper.s\\(x\\)", names(y))
    lowerCols <- grep("lower.s\\(x\\)", names(y))
    fit_ylim <- c(min(y[,lowerCols]), max(y[,upperCols]))
    if(!is.null(inputData)) {
      input_ylim <- c(min(sapply(inputData, min)), max(max(sapply(inputData, max))))
    }
  }
    
  if(!is.null(inputData)) {
    numTracks <- numTracks + ncol(inputData)
  }

  if(numTracks > 5) {
    par(mfrow = c(numTracks%/%2, numTracks%/%2 + numTracks%%2), oma = c(0, 0, 2, 0))
    xlab <- "genomic position"
  }
  else {
    par(mfrow = c(numTracks, 1), oma = c(0, 0, 2, 0))
    xlab <- ""
  }
  ## set penalty against scientific numbers
  opt <- options("scipen" = 100000)

  ## plot raw counts
  if(!is.null(inputData)) {
    for (ii in 1:ncol(inputData)) {
      plot(x, as.vector(inputData[,ii]), type = "p", col = "#73737330", pch = 19, 
           ylim = input_ylim, xlab = xlab, ylab = names(inputData)[ii])
    }
  }

  ## plot fits and confidence interval
  for(ii in sCols) {
    track <- names(y)[ii]

    ## print 'genomic position' only on the last plot
    if(ii == sCols[length(sCols)]) {
      xlab <- "genomic position"
    }
    ## if no ylim provided make sure that confidence intervalls are included
    if(!scale) {
      fit_ylim <- c(min(y[[paste("lower", track, sep = ".")]]),
                max(y[[paste("upper", track, sep = ".")]]))
    }

    plot(x, y[,ii], type = "l", col = 'black', ylim = fit_ylim, xlab = xlab, 
         ylab = track, ...)
    lines(x, y[[paste("lower", track, sep = ".")]], lty = "dotted")
    lines(x, y[[paste("upper", track, sep = ".")]], lty = "dotted")
    abline(h = 0, col = "red")
  }
  mtext(title, outer = TRUE, cex = 1.5)
  ## unset penalty for scientific numbers
  options(opt)
}
