## Get number of plot columns
getNumCols <- function(factorGroups) {
    ## all plot within groups
    pcols <- sum(sapply(factorGroups, function(y) {
      length(y)*(length(y) - 1)/2
    }))
    ## plus all plots between groups
    pcols <- pcols + length(factorGroups)*(length(factorGroups) - 1)/2
    return(pcols)
}

## plotting routine
plotQC <- function(object, factorGroups, plotPath) {
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
                    heatscatter(object[[ii]][,elem[jj]], object[[ii]][,elem[kk]], log = "xy",
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
                    heatscatter(x, y, log = "xy", main = paste("All", group1, "vs.", group2, names(object)[ii]),
                                xlab = group1, ylab = group2)
                    abline(0,1)
                }
            }
        }
    }
    dev.off()
}

## The quality check function for GAMDataSet
qcGenoGAMDataSet <- function(ggd, factorGroups = list()) {

    if(!requireNamespace("LSD", quietly = TRUE)) {
            stop("The 'LSD' package is required to plot proper heat-scatterplots. Please install it from CRAN", call. = FALSE)
    }
    
    qcPath <- file.path(getwd(), "qc")
    if(!dir.exists(qcPath)) {
        dir.create(qcPath)
    }
    hscatterFile <- "qc_normalization.png"

    countMat <- sum(ggd)[[1]]
    sf <- sizeFactors(ggd)
    countNorm <- t(t(countMat)/exp(sf))
    if(length(factorGroups) == 0) {
        factorGroups <- list(names(sf))
    }
    plotQC(list(normalized = countNorm, raw = countMat), factorGroups,
           plotPath = file.path(qcPath, hscatterFile))
}

## Find largest common substring
commonSubstring <- function(x) {
    res <- Biostrings::substr(x[1], start = 1, stop = lcprefix(x[1], x[length(x)]))
    if(res == "") {
        res <- paste(x, collapse = "_")
    }
    return(res)
}

## guess factorGroups from design matrix
guessFactorGroups <- function(design) {
    designCols <- sapply(design, is.integer)
    groups <- unique(design[,designCols, drop = FALSE])
    factorGroups <- lapply(1:nrow(groups), function(y) {
        commonNames <- as.list(design == groups[y,,drop = FALSE])[[1]]
        return(rownames(design)[commonNames])
    })
    return(factorGroups)
}
