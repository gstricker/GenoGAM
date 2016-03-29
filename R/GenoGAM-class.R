#############################
## GenoGAM class 
##################

#' @include GenoGAMSettings-class.R
#' @include GenomicTiles-class.R
NULL

#' GenoGAM class
#'
#' This class is designed to represent the model object containing the estimate parameters,
#' arguments and finals fits of the model on a basepair level.
#' 
#' @slot design A mgcv-type formula object.
#' @slot fits A data.frame of the fits, the standard error and the first and
#' second derivative of the fits for each experiment.
#' @slot positions A GPos object of the positions and seqnames corresponding
#' to the rows in the 'fits' slot.
#' @slot smooths A data.frame of knot positions and base function coefficients,
#' in order to reproduce the splines and compute derivatives. 
#' @slot experimentDesign The design matrix according to which the fitting
#' was performed.
#' @slot fitparams Global parameters 'lambda', 'theta', 'Coefficient of Variation' and
#' the 'penalty order' used to compute the model.
#' @slot family The distribution family.
#' @slot cvparams Parameters used for cross validation.
#' @slot settings The global and local settings that were used to compute the model.
#' @slot tileSettings A list of settings used to compute tiles.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @exportClass GenoGAM
setClass("GenoGAM",
         slots = list(design = "formula", fits = "data.frame",
             positions = "GPos", smooths = "data.frame",
             experimentDesign = "matrix", fitparams = "numeric",
             family = "ANY", cvparams = "numeric",
             settings = "GenoGAMSettings", tileSettings = "list"),
         prototype = prototype(design = ~ 1, fits = data.frame(),
             positions = GPos(), smooths = data.frame(),
             experimentDesign = matrix(), fitparams = numeric(), family = mgcv::nb(),
             cvparams = numeric(), settings = GenoGAMSettings(),
             tileSettings = list()))

## Validity
## ========

#' Validating the correct type
.validateDesignType <- function(object) {
    if(class(slot(object, "design")) != "formula") {
        return("'design' must be a formula object")
    }
    NULL
}

.validateFitsType <- function(object) {
    if(class(slot(object, "fits")) != "data.frame") {
        return("'fits' must be a data.frame object")
    }
    NULL
}

.validatePositionsType <- function(object) {
    if(class(slot(object, "positions")) != "GPos") {
        return("'positions' must be a GPos object")
    }
    NULL
}

.validateSmoothsType <- function(object) {
    if(class(slot(object, "smooths")) != "data.frame") {
        return("'smooths' must be a data.frame object")
    }
    NULL
}

.validateExpDesignType <- function(object) {
    if(class(slot(object, "experimentDesign")) != "matrix") {
        return("'experimentDesign' must be a matrix object")
    }
    NULL
}

.validateFitParamsType <- function(object) {
    if(class(slot(object, "fitparams")) != "numeric") {
        return("'fitparams' must be a numeric object")
    }
    NULL
}

.validateCVParamsType <- function(object) {
    if(class(slot(object, "cvparams")) != "numeric") {
        return("'cvparams' must be a numeric object")
    }
    NULL
}

.validateSettingsType <- function(object) {
    if(class(slot(object, "settings")) != "GenoGAMSettings") {
        return("'settings' must be a GenoGAMSettings object")
    }
    NULL
}

.validateTileSettingsType <- function(object) {
    if(class(slot(object, "tileSettings")) != "list") {
        return("'tileSettings' must be a list object")
    }
    NULL
}

## general validate function
.validateGenoGAM <- function(object) {
    c(.validateDesignType(object), .validateFitsType(object),
      .validatePositionsType(object), .validateSmoothsType(object),
      .validateExpDesignType(object), .validateFitParamsType(object),
      .validateCVParamsType(object), .validateSettingsType(object),
      .validateTileSettingsType(object))
}

setValidity2("GenoGAM", .validateGenoGAM)

## Constructor
## ========

#' GenoGAM constructor
#'
#' The GenoGAM constructor, not designed to be actually used, by the user.
#'
#' For arguments see slots of the class.
#' @return An object of the type GenoGAM.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.GenoGAM <- function(...) {
    return(new("GenoGAM", ...))
}

## Cosmetics
## ==========

#' The actual show function
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.showGenoGAM <- function(gg) {
    fitparams <- slot(gg, "fitparams")
    cvparams <- slot(gg, "cvparams")
    settings <- slot(gg, "tileSettings")
    
    show(slot(gg, "family"))
    cat("Formula:\n")
    show(slot(gg, "design"))
    cat("\n")
    
    cat("Experiment Design:\n")
    show(slot(gg, "experimentDesign"))
    cat("\n")
    cat("Global Estimates:\n")
    cat("  Lambda:", fitparams["lambda"], "\n")
    cat("  Theta:", fitparams["theta"], "\n")
    cat("  Coefficient of Variation:", fitparams["CoV"], "\n")

    cat("\n")
    cat("Cross Validation:")
    if(is.logical(cvparams["cv"]) & cvparams["cv"]) cat(" Performed\n")
    else cat(" Not performed\n")
    cat("  K-folds:", cvparams["kfolds"], "\n")
    cat("  Number of tiles:", cvparams["ncv"], "\n")
    cat("  Interval size:", cvparams["size"], "\n")

    cat("\n")
    cat("Tile settings:\n")
    cat("  chunk size:", settings$chunkSize, "\n")
    cat("  tile size:", settings$tileSize, "\n")
    cat("  overhang size:", settings$overhangSize, "\n")
    cat("  number of tiles:", settings$numTiles, "\n")
}

setMethod("show", "GenoGAM", function(object) {
    .showGenoGAM(object)
})

setMethod("summary", "GenoGAM", function(object) {
    fits <- data.table::data.table(cbind(seqnames(slot(object, "positions")),
                  GenomicRanges::pos(slot(object, "positions")),
                  slot(object, "fits")))
    data.table::setnames(fits, names(fits)[1:2], c("seqnames", "x"))
    .showGenoGAM(object)
    cat("\n")
    cat("Fitted values:\n")
    cat("\n")
    show(fits)
})

#' Make an example /code{GenoGAM}
#'
#' @return A /code{GenoGAM} object
#' @examples
#' test <- makeTestGenoGAM()
#' @export
makeTestGenoGAM <- function() {
    gg <- .GenoGAM()
    gp <- GenomicRanges::GPos(GenomicRanges::GRanges("chrI", IRanges::IRanges(1, 100)))
    fits <- data.frame("s(x)" = runif(100), "s(x):type" = runif(100))
    slot(gg, "positions") <- gp
    slot(gg, "fits") <- fits
    return(gg)
}

#' @rdname getFits
#' @export
setGeneric("getFits", function(object) standardGeneric("getFits"))

#' Accessor to \code{fits} slot
#'
#' The \code{fits} slot contains the fitted values of the model
#'
#' @param object A \code{GenomicTiles} object
#' @return A data.frame of the fits
#' @examples
#' gg <- makeTestGenoGAM()
#' fits <- getFits(gg)
#' @rdname getFits
#' @export
setMethod("getFits", "GenoGAM", function(object) object@fits)

.subsetByRanges <- function(gg, ranges) {
    pos <- slot(gg, "positions")
    ov <- findOverlaps(pos, ranges)
    indx <- queryHits(ov)
    slot(gg, "positions") <- pos[indx,]
    slot(gg, "fits") <- slot(gg, "fits")[indx,]
    return(gg)
}

.subsetByPosition <- function(gg, ...) {
    pos <- slot(gg, "positions")
    ov <- findOverlaps(pos, subset(pos, ...))
    indx <- queryHits(ov)
    slot(gg, "positions") <- pos[indx,]
    slot(gg, "fits") <- slot(gg, "fits")[indx,]
    return(gg)
}
#' Subset method for \code{GenoGAM}
#'
#' Subsetting the \code{GenoGAM} by a logical statement
#'
#' @param x A \code{GenoGAM} object.
#' @param ... Further arguments. Mostly a logical statement.
#' Note that the columnnames for chromosomes and positions
#' are: seqnames and pos.
#' @return A subsetted \code{GenoGAM} object.
#' @examples
#' gg <- makeTestGenoGAM()
#' subset(gg, pos <= 40)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setMethod("subset", "GenoGAM", function(x, ...) {
    .subsetByPosition(x, ...)
})

#' Subset by overlaps method for \code{GenoGAM}
#'
#' Subsetting the \code{GenoGAM} by a \code{GRanges} object
#'
#' @param query A \code{GenoGAM} object.
#' @param subject A \code{GRanges} object
#' @param ... Additional parameters
#' @return A subsetted \code{GenoGAM} object.
#' @examples
#' gg <- makeTestGenoGAM()
#' gr <- GRanges("chrI", IRanges(1,40))
#' subsetByOverlaps(gg, gr)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setMethod("subsetByOverlaps", "GenoGAM", function(query, subject) {
    .subsetByRanges(query, subject)
})

#' View the dataset
#'
#' Cbinding the columns all together and coercing to data.frame
#'
#' @param object A \code{GenoGAM} object
#' @param ranges A \code{GRanges} object. Makes it possible to
#' select regions by \code{GRanges}. Either ranges or seqnames, start and
#' end must be supplied
#' @param seqnames A chromosomes name. Either ranges or seqnames, start and
#' end must be supplied
#' @param start A start site. Either ranges or seqnames, start and
#' end must be supplied
#' @param end An end site. Either ranges or seqnames, start and
#' end must be supplied
#' @return A data.frame of the selected data.
#' @examples
#' gg <- makeTestGenoGAM()
#' gr <- GRanges("chrI", IRanges(1,40))
#' head(view(gg, gr))
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname GenoGAM-view
#' @export
setMethod("view", "GenoGAM", function(object, ranges = NULL, seqnames = NULL,
                                      start = NULL, end = NULL) {
    if((is.null(seqnames) | is.null(start) | is.null(end)) & is.null(ranges)) {
        temp <- object
    }
    else {
        if(!is.null(ranges)) {
            temp <- .subsetByRanges(object, ranges)
        }
        else {
            temp <- .subsetByPosition(object, seqnames == seqnames & pos >= start & pos <= end)
        }
    }
    res <- cbind(slot(temp, "positions"), slot(temp, "fits"))
    return(res)
})

.pvals <- function(gg, log.p = FALSE) {
    cols <- names(gg@fits)
    colindx <- which(unlist(regexec("se\\.", names(gg@fits))) < 0)

    res <- data.frame(matrix(NA, nrow(gg@fits), length(colindx)))
    for(ii in 1:length(colindx)) {
        colname <- cols[colindx[ii]]
        secolname <- paste("se", colname, sep = ".")
        background <- (regexec(":", colname)[[1]] < 0)
        if(background) {
            intercept <- median(gg@fits[,colname], na.rm = TRUE)
            fitMean <- abs(gg@fits[,colname] - intercept)
        }
        else fitMean <- abs(gg@fits[,colname])
        res[,ii] <- 2*pnorm(0, mean = fitMean, sd = gg@fits[, secolname], log.p = log.p)
    }
    names(res) <- paste("pvalue", cols[colindx], sep = ".")
    return(res)
}

.result <- function(gg, log.p = FALSE) {
    if(nrow(getFits(gg)) > 0) {
        pvals <- .pvals(gg, log.p)
        slot(gg, "fits") <- cbind(slot(gg, "fits"), pvals)
    }
    return(gg)
}

#' Compute significance.
#'
#' Based on the model fits this functions computes pointwise pvalues.
#'
#' @param gg A fitted GenoGAM object.
#' @param log.p Should pvalues be returned in log scale?
#' @return A GenoGAm object which fits has been updated by the pvalue columns.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
computeSignificance <- function(gg, log.p = FALSE) {
    .result(gg, log.p = log.p)
}


## #' Prediction from fitted GenoGAM
## #'
## #' Takes a fitted ‘GenoGAM’ object produced by ‘genogam()’ and returns
## #' predictions given a new set of values for the model covariates or
## #' the original values used for the model fit. Predictions can be
## #' accompanied by standard errors, based on the posterior distribution
## #' of the model coefficients.
## #'
## #' @param object A GenoGAM object.
## #' @param newdata A GenoGAMDataSet object parallel to the original data set.
## #' @param type Either 'terms', 'link' or 'response'. The first returns the
## #' fits for each spline function. The last two return the fit for the response
## #' variable transformed by the link function ('link') or not ('response').
## #' @param se.fit Should the standard errors be returned?
## #' @param terms Which terms should be returned?
## #' @param exclude Which terms should be excluded?
## #' @return A list of fits and optinally of standard errors.
## #' @author Georg Stricker \email{georg.stricker@@in.tum.de}
## #' @export
## predict.GenoGAM <- function(object, newdata, type = c("terms", "response", "link"),
##                             se.fit = FALSE,
##                             terms = NULL, exclude = NULL) {
##     type <- match.arg(type)
##     if(missing(newdata)) x <- slot(object, "positions")
##     else x <- rowRanges(newdata)
##     rows <- match(pos(x), pos(slot(object, "positions")))

##     secols <- which(unlist(regexec("se.fit", names(object@fits))) >= 1)

##     if(is.null(terms)) {
##         cols <- which(unlist(regexec("se.fit", names(object@fits))) != 1)
##     }
##     else {
##         cols <- sapply(terms, function(y) {
##             which(unlist(regexec(y, names(object@fits))) >= 1)
##         })
##         if(se.fit) secols <- cols[cols %in% secols]
##         cols <- cols[!(cols %in% secols)]
##     }

##     if(!is.null(exclude)) {
##       excols <- sapply(exclude, function(y) {
##           which(unlist(regexec(y, names(object@fits))) >= 1)
##       })
##       cols <- cols[!(cols %in% excols)]
##       secols <- secols[!(secols %in% excols)]
##     }
    
##     if(type == "response" | type == "link") {
##         cols <- which(unlist(regexec("se.fit", names(object@fits))) != 1)
##         fits <- rowSums(slot(object, "fits")[rows, cols])
##         if(se.fit) sefits <- rowSums(slot(object, "fits")[rows, secols])
##         if(type == "response") {
##             fits <- slot(object, "family")$linkinv(fits)
##             if(se.fit) sefits <- slot(object, "family")$linkinv(sefits)
##         }
##     }
##     else {
##         fits <- slot(object, "fits")[rows, cols]
##         if(se.fit) sefits <- slot(object, "fits")[rows, secols]
##     }
##     if(se.fit) return(list(fit = fits, se.fit = sefits))
##     return(list(fit = fits))
## }
