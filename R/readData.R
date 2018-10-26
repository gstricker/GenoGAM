### ====================
### read data functions
### ====================

## higher level read functions
## ===========================

#' Read Data function
#'
#' This is the core function to read and parse raw data from a config file.
#' At the moment only the BAM format is supported. It is not intended to be
#' used by the user directly, as it is called internally by the GenoGAMDataSet
#' constructor. However it is exported if people wish to separately assemble
#' their data and construct the GenoGAMDataSet from SummarizedExperiment
#' afterwards. It also offers the possibility to use the HDF5 backend.
#'
#' @param config A data.frame containing the experiment design of the model
#' to be computed with the first three columns fixed. See the 'experimentDesign'
#' parameter in \code{\link{GenoGAMDataSet}} or details here.
#' @param hdf5 Should the data be stored on HDD in HDF5 format? By default this
#' is disabled, as the Rle representation of count data already provides a
#' decent compression of the data. However in case of large organisms, a complex
#' experiment design or just limited memory, this might further decrease the
#' memory footprint.
#' @param split If TRUE the data will be stored as a list of DataFrames by
#' chromosome instead of one big DataFrame. This is only necessary if organisms
#' with a genome size bigger than 2^31 (approx. 2.14Gbp) are analyzed,
#' in which case Rs lack of long integers prevents having a well compressed Rle
#' of sufficient size.
#' @param settings A GenoGAMSettings object. Not needed by default, but might
#' be of use if only specific regions should be read in.
#' See \code{\link{GenoGAMSettings}}.
#' @param ... Further parameters that can be passed to low-level functions.
#' Mostly to pass arguments to custom process functions. In case the default
#' process functions are used, i.e. the default settings paramenter,
#' the most interesting parameters might be fragment length estimator method
#' from ?chipseq::estimate.mean.fraglen for single-end data.
#' @return A DataFrame of counts for each sample and position.
#' Or if split = TRUE, a list of DataFrames by chromosomes
#' @details
#' The config data.frame contains the actual experiment design. It must
#' contain at least three columns with fixed names: 'ID', 'file' and 'paired'.
#'
#' The field 'ID' stores a unique identifier for each alignment file.
#' It is recommended to use short and easy to understand identifiers because
#' they are subsequently used for labelling data and plots.
#'
#' The field 'file' stores the complete path to the BAM file.
#'
#' The field 'paired', values TRUE for paired-end sequencing data, and FALSE for
#' single-end sequencing data.
#'
#' Other columns will be ignored by this function.
#' @examples
#' # Read data
#' 
#' ## Set config file
#' config <- system.file("extdata/Set1", "experimentDesign.txt", package = "GenoGAM")
#' config <- read.table(config, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#' for(ii in 1:nrow(config)) {
#'     absPath <- system.file("extdata/Set1/bam", config$file[ii], package = "GenoGAM")
#'     config$file[ii] <- absPath
#' }
#'
#' ## Read all data
#' df <- readData(config)
#' df
#'
#' ## Read data of a particular chromosome
#' settings <- GenoGAMSettings(chromosomeList = "chrI")
#' df <- readData(config, settings = settings)
#' df
#' 
#' ## Read data of particular range
#' region <- GenomicRanges::GRanges("chrI", IRanges(10000, 20000))
#' params <- Rsamtools::ScanBamParam(which = region)
#' settings <- GenoGAMSettings(bamParams = params)
#' df <- readData(config, settings = settings)
#' df
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
readData <- function(config, hdf5 = FALSE, split = FALSE,
                     settings = GenoGAMSettings(), ...) {

    futile.logger::flog.info("Reading in data")
    
    fixedInput <- paste0("Using the following parameters:\n",
                    "  Split: ", split, "\n",
                    "  HDF5: ", hdf5, "\n")
    futile.logger::flog.debug(fixedInput)
    futile.logger::flog.debug(show(list(...)))
    futile.logger::flog.debug(show(settings))

    args <- list(...)

    ## read settings parameters
    ## get chromosomeLengths and check again chromosomeList
    header <- Rsamtools::scanBamHeader(config$file[1])
    chromosomeLengths <- header[[1]]$targets
    chromosomeList <- slot(settings, "chromosomeList")
    if(!is.null(chromosomeList)) {
        chromosomeLengths <- chromosomeLengths[names(chromosomeLengths) %in% chromosomeList]
        if(length(chromosomeLengths) == 0) {
            futile.logger::flog.error("The data does not match the region specification in the bamParams settings.")
            return(S4Vectors::DataFrame())
        }
    }
    ## get other parameters
    center <- slot(settings, "center")
    processFun <- slot(settings, "processFunction")
    if(is.null(processFun)) {
        slot(settings, "processFunction") <- .processCountChunks
        processFun <- slot(settings, "processFunction")
    }
    params <- slot(settings, "bamParams")

    ## read data
    rawData <- lapply(1:nrow(config), function(ii) {
        if(!is.null(center)) {
            args <- c(args, center = center)
        }
        futile.logger::flog.info(paste("Reading in", config$ID[ii]))
        futile.logger::flog.debug(paste(config$ID[ii], "is located at", config$file[ii],
                                        "and is paired end =", config$paired[ii]))
        res <- do.call(.readRawData, c(list(path = config$file[ii], processFUN = processFun,
                                            asMates = config$paired[ii], params = params,
                                            chromosomeLengths = chromosomeLengths), args))
        return(res)
    })
    names(rawData) <- config$ID

    ans <- lapply(rawData, unlist)
    names(ans) <- names(rawData)
    ans <- S4Vectors::DataFrame(ans)

    futile.logger::flog.info("Finished reading in data")
    return(ans)
}

## The intermediate function calling the correct
## function to read in data based on the suffix and
## computing the coverage from the GRanges output.
## Returns a list of coverage Rles per chromosome
.readRawData <- function(path, ...) {

    elements <- strsplit(path,"\\.")[[1]]
    suffix <- elements[length(elements)]
    res <- NULL
    
    if (suffix == "bam") {
        futile.logger::flog.debug(paste(path, "identified as BAM file"))
        res <- .readBAM(path = path, ...)
    }
    ## for the time being throw error if no BAM file
    else {
        stop("The input file doesn't seem to be in BAM format. At the moment only BAM format is supported")
    }

    grl <- GenomicRanges::GRangesList(res)
    GenomeInfoDb::seqlevels(grl) <- GenomeInfoDb::seqlevelsInUse(grl)
    coverageRle <- IRanges::coverage(grl)
    GenomeInfoDb::seqlengths(coverageRle) <- GenomeInfoDb::seqlengths(grl)

    if(length(Rsamtools::bamWhich(list(...)$params)) != 0) {
        chroms <- seqlevels(coverageRle)
        regions <- Rsamtools::bamWhich(list(...)$params)
        for(chr in chroms) {
            coverageRle[[chr]] <- coverageRle[[chr]][regions[[chr]]]
        }
    }

    return(coverageRle)
}


## Low level read in functions
## ===========================

#' A function to read BAM files.
#'
#' The functions reads in BAM files and processes them according to a given
#' process function. It returns a list of GRanges, one element per chromosome.
#'
#' @param path A character object indicating the path to the BAM file.
#' @param indexFile A character object indicating the path to the BAM Index
#' file. By default it is assumed to be in the same directory as the BAM
#' file.
#' @param processFUN A function specifying how to process the raw data.
#' @param chromosomeLengths A named numeric vector of chromosome lengths.
#' @param params A 'ScanBamParam' object defining the parameters to be used.
#' If NULL 'what' is defined as the columns 'pos' and 'qwidth'. 'Which' is
#' set according to the chunkSize parameter, but at most covers the entire
#' chromosome.
#' @param asMates A logical value indicating if mates are present in the BAM
#' file (paired-end) or not (single-end mostly).
#' @param ... Further parameters that are passed to the process functions.
#' @return A GRanges list of intervals
#' @noRd
.readBAM <- function(path, indexFile = path, processFUN, asMates,
                     chromosomeLengths, params = NULL, ...){

    ## stop if package for shifting in case of single-end data is not installed
    if(!asMates) {
        if(!requireNamespace("chipseq", quietly = TRUE)) {
            stop("The 'chipseq' package is required to estimate the fragment
length of single-end data. Please install it from Bioconductor or use
asMates = TRUE for paired-end data", call. = FALSE)
        }
    }

    args <- list(...)
    
    if(is.null(params)) {
        params <- Rsamtools::ScanBamParam(what = c("pos", "qwidth"))
    }
    
    ## convert which params to GRanges
    if(length(Rsamtools::bamWhich(params)) != 0) {
        regions <- GenomicRanges::GRanges(Rsamtools::bamWhich(params))
        lengths <- chromosomeLengths[GenomeInfoDb::seqlevels(regions)]
        suppressWarnings(GenomeInfoDb::seqlengths(regions) <- lengths)
        regions <- GenomicRanges::trim(regions)
        regions <- GenomicRanges::reduce(regions)
    }
    else {
        starts <- rep(1, length(chromosomeLengths))
        ends <- chromosomeLengths
        regions <- GenomicRanges::GRanges(names(chromosomeLengths),
                                          IRanges::IRanges(starts, ends))
    }

    futile.logger::flog.debug("The following regions will be read in:")
    futile.logger::flog.debug(show(regions))
    
    .local <- function(chromName, chromosomeCoords, asMates, path, indexFile,
                          params, processFUN, args) {
        ## load package for SnowParam or BatchJobs backend
        suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))

        ## read data
        coords <- chromosomeCoords[GenomeInfoDb::seqnames(regions) == chromName,]
        Rsamtools::bamWhich(params) <- coords
        
        if (asMates) {
            reads <- GenomicAlignments::readGAlignmentPairs(path, index = indexFile,
                                                            param = params)
        }
        else {
            reads <- GenomicAlignments::readGAlignments(path, index = indexFile,
                                                        param = params)
        }
        
        if(length(reads) == 0L) return(GRanges())
      
        ## ensure sequence levels are correct
        GenomeInfoDb::seqlevels(reads, pruning.mode="coarse") <- GenomeInfoDb::seqlevelsInUse(reads)
        
        ans <- do.call(processFUN, c(list(reads), args))
        return(ans)
    }
    
    ## run BAM reading in parallel. 
    res <- BiocParallel::bplapply(names(chromosomeLengths), .local,
                                  chromosomeCoords = regions, asMates = asMates, 
                                  path = path, indexFile = indexFile, params = params,
                                  processFUN = processFUN, args = args)
    
    return(res)
}

## Process functions for chunks of data
## ====================================

#' Processing raw data and convert to GRanges of intervals
#'
#' @param chunk An object of type GAlignments or GAlignmentPairs
#' @param center Logical value indicating if the fragments should be centered
#' or not. If centered, the fragment is reduced to only one data point,
#' which is the center of the fragment. If center = FALSE, then the complete
#' fragment is taken to compute the coverage. In single-end cases the
#' fragment size is estimated by one of the three methods: coverage,
#' correlation or SISSR. See ?chipseq::estimate.mean.fraglen.
#' @return A GRanges object of intervals.
#' @noRd
.processCountChunks <- function(chunk, center, ...) {
    paired <- switch(class(chunk),
                     GAlignmentPairs = TRUE,
                     GAlignments = FALSE)
    coords <- GenomicRanges::GRanges(GenomeInfoDb::seqlevels(chunk), IRanges::IRanges(1, GenomeInfoDb::seqlengths(chunk)))

    ## correct for unspecified reads and chunks containing only those
    if(!paired) {
        chunk <- chunk[IRanges::width(chunk) <= 2*median(IRanges::width(chunk))]
        strands <- droplevels(GenomicRanges::strand(chunk))
        if(!all(c("+", "-") %in% levels(strands))) {
            return(NULL)
        }
    }

    if(center) {
        fragments <- .centerFragments(chunk, asMates = paired, ...)
    }
    else {
        fragments <- .countFragments(chunk, asMates = paired, ...)
    }
    validFragments <- fragments[IRanges::start(fragments) >= IRanges::start(coords) &
                                IRanges::end(fragments) <= IRanges::end(coords)]

    return(fragments)
}

#' A function to center fragments.
#'
#' This functions centers fragments with regard to the strand and returns a
#' GRanges object with centered ranges, that is ranges of width one.
#'
#' @param reads A GAlignment object as returned by the 'readGAlignments'
#' functions.
#' @param asMates A logical value indicating if mates are present in the BAM
#' file (paired-end) or not (single-end mostly).
#' @param shiftMethods The method to be used when estimating the fragment
#' length for single-end reads (see ?chipseq::estimate.mean.fraglen).
#' Other methods are 'SISSR' and 'correlation'.
#' @param ... Further parameters that can be passed on to
#' chipseq::estimate.mean.fraglen.
#' @return A GRanges object of centered fragments.
#' @noRd
.centerFragments <- function(reads, asMates, shiftMethod = c("coverage", "correlation", "SISSR"), ...) {
    if(is.null(reads)) return(NULL)

    plusStrand <- reads[GenomicRanges::strand(reads) == "+",]
    negStrand <- reads[GenomicRanges::strand(reads) == "-",]
    
    if (asMates) {
        plusMidpoints <- .getPairedCenters(plusStrand)
        negMidpoints <- .getPairedCenters(negStrand)
    }
    else {

        ## estimate fragment length for shift
        granges_reads <- GenomicRanges::granges(reads)
        fraglen <- chipseq::estimate.mean.fraglen(granges_reads, method =
                                                  match.arg(shiftMethod),
                                                  ...)
        plusMidpoints <- GenomicRanges::granges(plusStrand)
        IRanges::start(plusMidpoints) <- IRanges::end(plusMidpoints) <- IRanges::start(plusMidpoints) + fraglen/2
        
        negMidpoints <- GenomicRanges::granges(negStrand)
        IRanges::end(negMidpoints) <- IRanges::start(negMidpoints) <- IRanges::end(negMidpoints) - fraglen/2
    }

    midpoints <- sort(c(plusMidpoints,negMidpoints))
    return(midpoints)
}

#' A function to get the fragment centers
#'
#' @param chunk A GAligmnents or GAlignmentPairs object.
#' @return A GRanges object of centered fragments, that is of lenght one
#' @noRd
.getPairedCenters <- function(chunk) {
    if(length(chunk) == 0) return(GenomicRanges::GRanges())
    gr <- GenomicRanges::GRanges(chunk)

    fragmentSize <- IRanges::width(gr)
    midpoints <- (IRanges::start(gr) + IRanges::end(gr))/2
  
    IRanges::start(gr) <- IRanges::end(gr) <- midpoints

    ## filter for artefacts because of wrongly mapped reads
    maxAllowedFragSize <- 2*median(fragmentSize)
    falseFragments <- which(fragmentSize > maxAllowedFragSize)
    
    futile.logger::flog.debug(paste(length(falseFragments), "dropped due to maximum allowed fragment size of:",
                                    maxAllowedFragSize, "on the following GRange:"))
    futile.logger::flog.debug(show(gr))
    
    if (length(falseFragments) > 0) {
        gr <- gr[-falseFragments]
    }
    return(gr)
}

#' Extracting fragments given reads from a BAM file.
#'
#' @param reads A GAligments or GAlignmenPairs object as returned by the 
#' 'readGAlignments' functions.
#' @param asMates A logical value indicating if mates are present in the BAM
#' file (paired-end) or not (single-end mostly)
#' @param shiftMethod The method to be used when estimating the fragment
#' length for single-end reads (see ?chipseq::estimate.mean.fraglen). Other
#' methods are 'SISSR' and 'correlation'.
#' @param ... Further parameters that can be passed on to.
#' chipseq::estimate.mean.fraglen.
#' @return A GRanges object of full fragments
#' @noRd
.countFragments <- function(reads, asMates, shiftMethod =
                            c("coverage", "correlation", "SISSR"), ...) {
    if(is.null(reads)) return(NULL)
    
    plusStrand <- reads[GenomicRanges::strand(reads) == "+",]
    negStrand <- reads[GenomicRanges::strand(reads) == "-",]

    if (asMates) {
        plusFragments <- .getFragments(plusStrand)   
        negFragments <- .getFragments(negStrand)
    }
    else {
        granges_reads <- GenomicRanges::granges(reads)
        fraglen <- chipseq::estimate.mean.fraglen(granges_reads, method =
                                                  match.arg(shiftMethod),
                                                  ...)
        plusFragments <- GenomicRanges::granges(plusStrand)
        IRanges::end(plusFragments) <- IRanges::start(plusFragments) + fraglen
        
        negFragments <- GenomicRanges::granges(negStrand)
        IRanges::start(negFragments) <- IRanges::end(negFragments) - fraglen
    }

    allFragments <- sort(c(plusFragments,negFragments))
    return(allFragments)
}

#' A function to get the fragments
#'
#' @param chunk A GAligmnents or GAlignmentPairs object.
#' @return A GRanges object of fragments
#' @noRd
.getFragments <- function(chunk){
    if(length(chunk) == 0) return(GenomicRanges::GRanges())
    gr <- GenomicRanges::GRanges(chunk)

    maxAllowedFragSize <- 2*median(IRanges::width(gr))
    falseFragments <- which(IRanges::width(gr) > maxAllowedFragSize)

    futile.logger::flog.debug(paste(length(falseFragments), "dropped due to maximum allowed fragment size of:",
                                    maxAllowedFragSize, "on the following GRange:"))
    futile.logger::flog.debug(show(gr))
    
    if (length(falseFragments) > 0) {
        gr <- gr[-falseFragments]
    }
    return(gr)
}
