### ===========================================================
### Genomic Tiles class
### ===========================================================

#' @include GenoGAMSettings-class.R
#' @include GenomicTiles-class.R
NULL

#' GenoGAMDataSet
#'
#' This class is designed to represent the input for the GenoGAM model.
#' it extends the GenomicTiles class.
#' 
#' @slot settings The global and local settings that were used to compute the model.
#' @slot design The formula describing how to evaluate the data.
#' @slot sizeFactors The normalized values for each sample. 
#' @details For all other slots see SummarizedExperiment.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @exportClass GenoGAMDataSet
setClass("GenoGAMDataSet",
         contains = "GenomicTiles",
         slots = list(settings = "GenoGAMSettings",
             design = "formula",
             sizeFactors = "numeric"),
         prototype = list(settings = GenoGAMSettings(),
             design = ~ 1,
             sizeFactors = numeric()))

## Validity
## ========

.validateSettingsType <- function(object) {
    if(class(slot(object, "settings")) != "GenoGAMSettings") {
        return("'settings' must be a GenoGAMSettings object")
    }
    NULL
}

.validateDesignType <- function(object) {
    if(class(slot(object, "design")) != "formula") {
        return("'design' must be a formula object")
    }
    NULL
}

.validateSFType <- function(object) {
    if(class(slot(object, "sizeFactors")) != "numeric") {
        return("'sizeFactors' must be a numeric object")
    }
    NULL
}

## general validate function
.validateGenoGAMDataSet <- function(object) {
    c(.validateSettingsType(object), .validateDesignType(object),
      .validateSFType(object))
}

setValidity2("GenoGAMDataSet", .validateGenoGAMDataSet)

## Constructor
## ===========

#' GenoGAMDataSet constructor.
#'
#' This is the constructor function for GenoGAMDataSet. So far a GenoGAMDataSet
#' can be constructed from either an experiment design file or data.frame or
#' directly from a RangedSummarizedExperiment with a GPos object being the rowRanges.
#'
#' @param experimentDesign Either a character object specifying the path to a
#' delimited text file (the delimiter will be determined automatically),
#' or a data.frame specifying the experiment design. See details for the
#' structure of the experimentDesign.
#' @param chunkSize An integer specifying the size of one chunk in bp.
#' @param overhangSize An integer specifying the size of the overhang in bp.
#' As the overhang is taken to be symmetrical, only the overhang of one side
#' should be provided.
#' @param design A mgcv-like formula object. See details for its structure.
#' @param directory The directory from which to read the data. By default
#' the current working directory is taken.
#' @param settings A GenoGAMSettings object. This class is already present but 
#' not yet fully tested and therefore not accessible to the user. This
#' argument exists however in order to allow some workarounds if necessary.
#' See the vignette for a possible use.
#' @param ... Further parameters, mostly for arguments of custom processing
#' functions or to specify a different method for fragment size estimation.
#' See details for further information.
#' @return An object of class GenoGAMDataSet.
#' @details The experimentDesign file/data.frame must contain at least three
#' columns with fixed names: 'ID', 'file' and 'paired'.The field 'ID' stores
#' a unique identifier for each alignment file. It is recommended to use short
#' and easy to understand identifiers because they are subsequently used for
#' labelling data and plots. The field 'file' stores the BAM file name.
#' The field 'paired', values TRUE for paired-end sequencing data, and FALSE for
#' single-end sequencing data.  All other columns are stored in the colData
#' slot of the GenoGAMDataSet object. Note that all columns which will be used for
#' analysis must have at most two conditions, which are for now restricted
#' to 0 and 1. For example, if the IP data schould be corrected for input,
#' then the input will be 0 and IP will be 1, since we are interested in the
#' corrected IP. See examples.
#'
#' Design must be a mgcv-like formula. At the moment only the following is
#' possible: Either '~ 1' for a constant. ~ s(x) for a smooth fit over the
#' entire data. s(x, by = "myColumn"), where 'myColumn' is a column name
#' in the experimentDesign. This type of formula will then only fit the
#' samples annotated with 1 in this column.
#' Or ~ s(x) + s(x, by = "myColumn") + s(x, by = ...) + ...
#' The last formula lets you combine any number of columns, given they are
#' binary with 0 and 1. For example the formula for correcting IP for
#' input would look like this: ~ s(x) + s(x, by = "experiment"), where
#' 'experiment' is a column with 0s and 1s, with the ip samples annotated
#' with 1 and input samples with 0.
#''
#' In case of single-end data in might be usefull to specify a different
#' method for fragment size estimation. The argument 'shiftMethod' can be
#' supplied with the values 'coverage' (default), 'correlation' or 'SISSR'.
#' See ?chipseq::estimate.mean.fraglen for explanation.
#' @examples
#' \dontrun{
#' myConfig <- data.frame(ID = c("input","ip"),
#'                   file = c("myInput.bam", "myIP.bam"),
#'                   paired = c(FALSE, FALSE),
#'                   experiment = factor(c(0,1)),
#'                   stringsAsFactors = FALSE) 
#' myConfig2 <- data.frame(ID = c("wildtype1","wildtype2",
#'                               "mutant1", "mutant2"),
#'                   file = c("myWT1.bam", "myWT2.bam"
#'                            "myMutant1.bam", "myMutant2.bam"),
#'                   paired = c(FALSE, FALSE, FALSE, FALSE),
#'                   experiment = factor(c(0, 0, 1, 1)),
#'                   stringsAsFactors = FALSE)
#' 
#' gtiles <- GenoGAMDataSet(myConfig, chunkSize = 2000,
#' overhang = 250, design = ~ s(x) + s(x, by = "experiment")
#' gtiles <- GenoGAMDataSet(myConfig2, chunkSize = 2000,
#' overhang = 250, design = ~ s(x) + s(x, by = "experiment"))
#' }
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
GenoGAMDataSet <- function(experimentDesign, chunkSize, overhangSize, design,
                           directory = ".", settings = NULL, ...) {

    if(missing(experimentDesign)) {
        gt <- GenomicTiles()
        return(new("GenoGAMDataSet", gt))
    }

    if(class(experimentDesign) == "RangedSummarizedExperiment") {
        gt <- .GenoGAMDataSetFromSE(se = experimentDesign,
                                    chunkSize = chunkSize,
                                    overhangSize = overhangSize,
                                    design = design,
                                    settings = settings, ...)
    }
    else {
        gt <- .GenoGAMDataSetFromConfig(config = experimentDesign,
                                        chunkSize = chunkSize,
                                        overhangSize = overhangSize,
                                        design = design,
                                        directory = directory,
                                        settings = settings, ...)
    }
    
    return(gt)
}

#' Convert the config columns to the right type.
#'
#' @param config A data.frame with pre-specified columns.
#' @return The same data.frame with the columns of the right type.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.normalizeConfig <- function(config, directory) {
    
    if(class(config) == "character") {
        config <- fread(config, header = TRUE, data.table = FALSE)
    }

    config$ID <- as.factor(config$ID)
    config$file <- file.path(directory, as.character(config$file))
    config$paired <- as.logical(config$paired)
   
    return(config)
}

#' Construct GenomicTiles from a config file or a config data.frame
#'
#' See GenomicTiles in GenomicTiles-class.R for description.
#' @return A GenomicTiles object.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.GenoGAMDataSetFromConfig <- function(config, chunkSize, overhangSize,
                                    design, directory, settings, ...) {

    ## initialize some variables
    args <- list()
        
    ## normalize config object
    config <- .normalizeConfig(config, directory)

    ## get default settings
    if(is.null(settings)) settings <- GenoGAMSettings()
    center <- getDefaults(settings, "center")
    if(!is.null(center)) {
        settings <- setDefaults(settings,
                                processFunction = .processCountChunks)
    }

    ## get chromosomeLengths
    header <- Rsamtools::scanBamHeader(config$file[1])
    chroms <- header[[1]]$targets
    chromosomeLengths <- GenomicRanges::GRanges(names(chroms), IRanges(start = rep(1, length(chroms)),
                                                        end = chroms))
    GenomeInfoDb::seqlengths(chromosomeLengths) <- chroms
    chromosomeList <- getDefaults(settings, "chromosomeList")
    if(!is.null(chromosomeList)) {
        GenomeInfoDb::seqlevels(chromosomeLengths, force = TRUE) <- chromosomeList
        chroms <- chroms[names(chroms) %in% chromosomeList]
    }

    ## read in data
    futile.logger::flog.info("Reading in data.")
    rawData <- lapply(1:nrow(config), function(ii) {
        if(!is.null(center)) {
            args <- list(paired = config$paired[ii],
                         center = center)
        }
        unlist(do.call(.readData,
                      c(list(path = config$file[ii],
                             processFUN =
                             getDefaults(settings, "processFunction"),
                             chromosomeList = chromosomeList,
                             params = getDefaults(settings, "bamParams"),
                             asMates = config$paired[ii]),
                        args, list(...))), use.names = FALSE)
    })

    names(rawData) <- config$ID
    assays <- DataFrame(rawData)

    ## generate rowRanges
    bamParamsWhich <- Rsamtools::bamWhich(getDefaults(settings, "bamParams"))
    if(length(bamParamsWhich) != 0) {
        gp <- GenomicRanges::GPos(GenomicRanges::GRanges(bamParamsWhich))
        GenomeInfoDb::seqlengths(gp) <- chroms[GenomeInfoDb::seqlevels(GRanges(bamParamsWhich))]
    }
    else {
        gp <- GenomicRanges::GPos(chromosomeLengths)
        GenomeInfoDb::seqlengths(gp) <- chroms
    }
    
    gt <- GenomicTiles(assays, rowRanges = gp, chunkSize = chunkSize,
                       overhangSize = overhangSize)
    metadata(slot(gt, "index"))$check <- TRUE

    ## make colData
    colData <- DataFrame(config)[,-c(1:3), drop = FALSE]
    rownames(colData) <- config$ID
    colData(gt) <- colData
    ## initiate size factors
    sf <- rep(0, ncol(gt))
    names(sf) <- colnames(gt)

    settings <- setDefaults(settings, chromosomeList = GenomeInfoDb::seqlevels(tileSettings(gt)$chromosomes))

    gtd <- new("GenoGAMDataSet", gt, settings = settings,
               design = design, sizeFactors = sf)

    ## check if everything was set fine
    correct <- checkSettings(gtd)
    if(!correct) break
    
    futile.logger::flog.info("DONE")
    return(gtd)
}

.GenoGAMDataSetFromSE <- function(se, chunkSize, overhangSize,
                                    design, settings, ...) {
    if(is.null(settings)) settings <- GenoGAMSettings()

    gt <- GenomicTiles(se, chunkSize = chunkSize,
                       overhangSize = overhangSize)
    metadata(slot(gt, "index"))$check <- TRUE

    sf <- rep(0, ncol(gt))
    names(sf) <- colnames(gt)

    settings <- setDefaults(settings, chromosomeList = GenomeInfoDb::seqlevels(tileSettings(gt)$chromosomes))

    gtd <- new("GenoGAMDataSet", gt, settings = settings,
               design = design, sizeFactors = sf)

    ## check if everything was set fine
    correct <- checkSettings(gtd)
    if(!correct) break
    
    return(gtd)   
}

#' Make an example /code{GenoGAMDataSet}
#'
#' @return A /code{GenoGAMDataSet} object
#' @examples
#' test <- makeTestGenoGAMDataSet()
#' @export
makeTestGenoGAMDataSet <- function() {
    gp <- GenomicRanges::GPos(GenomicRanges::GRanges(c("chrI", "chrII"), IRanges(c(1,1), c(50,50))))
    df <- DataFrame(a = Rle(1:100), b = Rle(101:200))
    se <- SummarizedExperiment::SummarizedExperiment(list(df), rowRanges = gp)
    ggd <- GenoGAMDataSet(se, chunkSize = 15, overhangSize = 3,
                          design = ~s(x))
}

## Accessors
## =========

#' Access the \code{design} slot
#'
#' The \code{design} slot contains the \code{formula}
#' object which is used to fit the model
#'
#' @param object A GenoGAMDataSet object.
#' @param value A \code{formula} object
#' @return A \code{formula} object
#' @examples
#' ggd <- makeTestGenoGAMDataSet()
#' design(ggd)
#' design(ggd) <- ~1
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
#' @rdname design
setMethod("design", "GenoGAMDataSet", function(object) {
    slot(object, "design")
})

#' @export
#' @rdname design
setReplaceMethod("design", "GenoGAMDataSet", function(object, value) {
    slot(object, "design") <- value
    return(object)
})

#' Access the \code{sizeFactor} slot
#'
#' The \code{sizeFactor} slot contains the vector of
#' normalization values for each sample
#'
#' @param object A GenoGAMDataSet object.
#' @param value A named numeric vector
#' @return A named numeric vector
#' @examples
#' ggd <- makeTestGenoGAMDataSet()
#' sizeFactors(ggd)
#' sizeFactors(ggd) <- c(a = 5, b = 1/5)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
#' @rdname sizeFactors
setMethod("sizeFactors", "GenoGAMDataSet", function(object) {
    slot(object, "sizeFactors")
})

#' @export
#' @rdname sizeFactors
setReplaceMethod("sizeFactors", "GenoGAMDataSet", function(object, value) {
    slot(object, "sizeFactors") <- value
    return(object)
})

## Coercion
## ========

# converts the GenomicTiles object to a DataFrame
.GenoGAMDataSetToDataFrame <- function(from) {
    df <- as.data.frame(rowRanges(from))
    numDataFrames <- length(assays(from))
    gtDim <- dim(from)
    
    res <- unlist(assays(from))   
    gtdf <- DataFrame(cbind(df, res))

    metadata(gtdf) <- list(sizeFactors = sizeFactors(from))
        
    return(gtdf)
}

#' GenoGAMDataSet to DataFrame
#' 
#' @name GenoGAMDataSetToDataFrame
#' @family GenoGAMDataSet
setAs(from = "GenoGAMDataSet", to = "DataFrame", def = function(from) {
    .GenoGAMDataSetToDataFrame(from)
})

setMethod("as.data.frame", "GenoGAMDataSet", function(x) {
    ## Bug in as(x, "data.frame"), replace with
    as.data.frame(as(x, "DataFrame"))
})

## Subsetting
## ==========
#' Subset method for \code{GenoGAMDataSet}
#'
#' Subsetting the \code{GenoGAMDataSet} by a logical statement
#'
#' @param x A \code{GenoGAMDataSet} object.
#' @param ... Further arguments. Mostly a logical statement.
#' Note that the columnnames for chromosomes and positions
#' are: seqnames and pos.
#' @return A subsetted \code{GenomicTiles} object.
#' @examples
#' ggd <- makeTestGenoGAMDataSet()
#' res <- subset(ggd, seqnames == "chrI" & pos <= 50)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setMethod("subset", "GenoGAMDataSet", function(x, ...) {
    settings <- getSettings(x)
    design <- design(x)
    sf <- sizeFactors(x)
    gt <- .subsetGenomicTiles(x, ...)
    settings <- setDefaults(settings, chromosomeList = GenomeInfoDb::seqlevels(gt))

    gtd <- new("GenoGAMDataSet", gt, settings = settings,
               design = design, sizeFactors = sf)
    return(gtd)
})

#' Subset by overlaps method for \code{GenoGAMDataSet}
#'
#' Subsetting the \code{GenoGAMDataSet} by a \code{GRanges} object
#'
#' @param query A \code{GenoGAMDataSet} object.
#' @param subject A \code{GRanges} object
#' @param maxgap,minoverlap Intervals with a separation of \code{maxgap} or
#' less and a minimum of \code{minoverlap} overlapping positions, allowing for
#' \code{maxgap}, are considered to be overlapping.  \code{maxgap} should
#' be a scalar, non-negative, integer. \code{minoverlap} should be a scalar,
#' positive integer.
#' @param type By default, any overlap is accepted. By specifying the \code{type}
#' parameter, one can select for specific types of overlap. The types correspond
#' to operations in Allen's Interval Algebra (see references). If \code{type}
#' is \code{start} or \code{end}, the intervals are required to have matching
#' starts or ends, respectively. While this operation seems trivial, the naive
#' implementation using \code{outer} would be much less efficient. Specifying
#' \code{equal} as the type returns the intersection of the \code{start} and
#' \code{end} matches. If \code{type} is \code{within}, the query interval must
#' be wholly contained within the subject interval. Note that all matches must
#' additionally satisfy the \code{minoverlap} constraint described above.
#'
#' The \code{maxgap} parameter has special meaning with the special
#' overlap types. For \code{start}, \code{end}, and \code{equal}, it specifies
#' the maximum difference in the starts, ends or both, respectively. For
#' \code{within}, it is the maximum amount by which the query may be wider
#' than the subject.
#' @param ... Additional parameters
#' @return A subsetted \code{GenoGAMDataSet} object.
#' @examples
#' ggd <- makeTestGenoGAMDataSet()
#' gr <- GRanges("chrI", IRanges(1,50))
#' res <- subsetByOverlaps(ggd, gr)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setMethod("subsetByOverlaps", c("GenoGAMDataSet", "GRanges"),
          function(query, subject, maxgap=0L, minoverlap=1L,
                      type=c("any", "start", "end", "within", "equal"),...) {
              settings <- getSettings(query)
              design <- design(query)
              sf <- sizeFactors(query)
              subgt <- .subsetByOverlaps(query, subject, maxgap=maxgap,
                                         minoverlap=minoverlap,
                                         type=type,...)
              settings <- setDefaults(settings, chromosomeList = GenomeInfoDb::seqlevels(subgt))
              
              gtd <- new("GenoGAMDataSet", subgt, settings = settings,
                         design = design, sizeFactors = sf)
              return(gtd)
          })

#' Subsetting by GRanges
#'
#' Providing subsetting by GRanges through the single-bracket operator
#'
#' @param x A \code{GenoGAMDataSet} object
#' @param i A \code{GRanges} object
#' @return A subsetted \code{GenoGAMDataSet} object
#' @rdname GenoGAMDataSet-brackets
setMethod("[", c("GenoGAMDataSet", "GRanges"), function(x, i) {
    settings <- getSettings(x)
    design <- design(x)
    sf <- sizeFactors(x)
    subgt <- .exactSubsetByOverlaps(x, i)
    settings <- setDefaults(settings, chromosomeList = GenomeInfoDb::seqlevels(subgt))
              
    gtd <- new("GenoGAMDataSet", subgt, settings = settings,
               design = design, sizeFactors = sf)
    return(gtd)
})
    

## Silent methods
## ==============
setGeneric("getSettings", function(object) standardGeneric("getSettings"))

setMethod("getSettings", "GenoGAMDataSet", function(object) {
    slot(object, "settings")
})
