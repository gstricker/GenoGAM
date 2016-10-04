## Functions for quality checks for different objects

## Get number of plot columns dependent on the factor groups
#' @param factorGroups A list representing the groups. Each element is a group containing
#' the names of the samples belonging to it
#' @return Number of columns needed for multi-plotting
#' @noRd
getNumCols <- function(factorGroups) {
    ## all plot within groups
    pcols <- sum(sapply(factorGroups, function(y) {
      length(y)*(length(y) - 1)/2
    }))
    ## plus all plots between groups
    pcols <- pcols + length(factorGroups)*(length(factorGroups) - 1)/2
    return(pcols)
}

## plotting routine for GenoGAMDataSet
#' @param object A \code{GenoGAMDataSet} object
#' @param factorGroups A list representing the groups. Each element is a group containing
#' the names of the samples belonging to it
#' @param plotPath The file to plot to
#' @return NULL. A plot is saved to given directory
#' @noRd
plotQC_GenoGAMDataSet <- function(object, factorGroups, plotPath) {
    pcols <- getNumCols(factorGroups)
    prows <- length(object)
    png(plotPath, height = 1000, width = 1600)
    par(mfrow = c(prows, pcols))
    par(mar=c(5,4,4,1))
    for(ii in 1:length(object)) {
        counter <- 0
        while(counter < length(factorGroups)) {
            elem <- factorGroups[[counter + 1]]
            len <- length(elem)
            if(len <= 1) {
                counter <- counter + 1
                next
            }
            for (jj in 1:(len - 1)) {
                for(kk in (jj + 1) : len) {
                    LSD::heatscatter(object[[ii]][,elem[jj]], object[[ii]][,elem[kk]], log = "xy",
                                main = paste(elem[jj], "vs.", elem[kk], names(object)[ii]),
                                xlab = elem[jj], ylab = elem[kk])
                    abline(0,1)
                }
            }
            counter <- counter + 1
        }
        if(length(factorGroups) > 1) {
            lenGroups <- length(factorGroups)
            for(ll in 1:(lenGroups -1)) {
                for(mm in (ll + 1): lenGroups) {
                    group1 <- commonSubstring(factorGroups[[ll]])
                    group2 <- commonSubstring(factorGroups[[mm]])
                    if(length(factorGroups[[ll]]) > 1) {
                        x <- rowMeans(object[[ii]][,factorGroups[[ll]]])
                    }
                    else {
                        x <- object[[ii]][,factorGroups[[ll]]]
                    }
                    if(length(factorGroups[[mm]]) > 1) {
                        y <- rowMeans(object[[ii]][,factorGroups[[mm]]])
                    }
                    else {
                        y <- object[[ii]][,factorGroups[[mm]]]
                    }
                    LSD::heatscatter(x, y, log = "xy", main = paste("All", group1, "vs.", group2, names(object)[ii]),
                                xlab = group1, ylab = group2)
                    abline(0,1)
                }
            }
        }
    }
    suppressMessages(dev.off())
}

## plot histogram of counts
#' @param object A \code{GenoGAMDataSet}
#' @param plotPath The file to plot to
#' @return Nothing. Saves plot in plotPath
#' @noRd
plotQC_hist <- function(object, plotPath) {
    numVars <- ncol(object)
    totalSums <- colSums(object)
        
    png(plotPath, height = 800, width = 1200)
    par(mfrow = c(numVars, 1))
    for(ii in 1:numVars) {
        lab <- colnames(object)[ii]
        prettyTotal <- format(totalSums[lab], big.mark = ",", scientific = FALSE)
        hist(object[,ii], breaks = 100, xlab = lab,
             main = paste0("Total counts in ", lab, ": ", prettyTotal))
        abline(v = median(object[,ii]), col = "red", lwd = 2)
        axis(side = 1, at = median(object[,ii]), col.axis = "red")
    }
    suppressMessages(dev.off())
}

## The quality check function for GenoGAMDataSet
#' @param ggd A \code{GenoGAMDataSet} object
#' @param factorGroups A list representing the groups. Each element is a group containing
#' the names of the samples belonging to it
#' @param name A name specification for the file. Default is 'qc'
#' @return NULL. A plot is saved to given directory
#' @noRd
qcGenoGAMDataSet <- function(ggd, factorGroups = list(), name = 'qc') {

    if(!requireNamespace("LSD", quietly = TRUE)) {
            stop("The 'LSD' package is required to plot proper heat-scatterplots. Please install it from CRAN", call. = FALSE)
    }
    
    qcPath <- file.path(getwd(), "qc")
    if(!dir.exists(qcPath)) {
        dir.create(qcPath)
    }

    ## generate file name
    timestamp <- as.character(Sys.time())
    timestamp <- gsub("-", "", timestamp)
    timestamp <- gsub(":", "", timestamp)
    timestamp <- gsub(" ", "_", timestamp)
    hscatterFile <- paste0(name, "_normalization_", timestamp, ".png")
    densityFile <- paste0(name, "_counts_", timestamp, ".png")

    ## prepare count matrix
    countMat <- sum(ggd)[[1]]
    sf <- sizeFactors(ggd)
    countNorm <- t(t(countMat)/exp(sf))
    if(length(factorGroups) == 0) {
        factorGroups <- list(names(sf))
    }

    plotQC_hist(countMat, plotPath = file.path(qcPath, densityFile))
    
    plotQC_GenoGAMDataSet(list(normalized = countNorm, raw = countMat), factorGroups,
                          plotPath = file.path(qcPath, hscatterFile))
    message(paste0("Plots saved in\n", file.path(qcPath, densityFile), "\n",
                  file.path(qcPath, hscatterFile)))
}

#' A function to quality check the data
#'
#' This function checks some data attributes in the given class.
#' Check details for more information.
#'
#' @param object Any object for which this methods is implemented
#' @param ... further parameters. See details.
#' @return Based on the object provided, see details.
#' @details So far this method is only implemented for the class \code{GenoGAMDataSet}.
#' In this case some general metrics are printed and some plots are stored in the folder
#' "qc", which will be created in the working directory.
#'
#' Additional parameters:
#' factorGroups (for GenoGAMDataSet), which is used to specify factor groups for normalization plots.
#' By default the groups will be identified automatically. See ?computeSizeFactors for parameter description.
#'
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
qualityCheck <- function(object, ...) {
    if(class(object) == "GenoGAMDataSet") {
       res <- qcGenoGAMDataSet(object, ...)
    }
    return(res)
}

## Find largest common substring in names, based on common naming pattern for replicates
## this does not work for all situations.
#' @param x A character vector
#' @return The common largest substring
#' @noRd
commonSubstring <- function(x) {
    res <- Biostrings::substr(x[1], start = 1, stop = lcprefix(x[1], x[length(x)]))
    if(res == "") {
        res <- paste(x, collapse = "_")
    }
    return(res)
}

## guess factorGroups from design matrix
#' @param design The design matrix as provided by colData(GenoGAMDataSet)
#' @return A list of factor groups
#' @noRd
guessFactorGroups <- function(design) {
    designCols <- sapply(design, is.integer)
    groups <- unique(design[,designCols, drop = FALSE])
    factorGroups <- lapply(1:nrow(groups), function(y) {
        commonNames <- as.list(design == groups[y,,drop = FALSE])[[1]]
        return(rownames(design)[commonNames])
    })
    return(factorGroups)
}
