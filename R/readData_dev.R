### ====================
### read data functions
### ====================

## Reading in files
## ================

#' A function to read BAM files.
#'
#' The functions reads in BAM files and processes them according to a given
#' process function. It returns an RleList with the chromosome RLEs as
#' elements.
#'
#' @param path A character object indicating the path to the BAM file.
#' @param indexFile A character object indicating the path to the BAM Index
#' file. By default it is assumed to be in the same directory as the BAM
#' file.
#' @param processFUN A function specifying how to process the raw data.
#' @param chromosomeList A character vector of chromosome names. If NULL the
#' list is taken from the header of the BAM file.
#' @param params A 'ScanBamParam' object defining the parameters to be used.
#' If NULL 'what' is defined as the columns 'pos' and 'qwidth'. 'Which' is
#' set according to the chunkSize parameter, but at most covers the entire
#' chromosome.
#' @param asMates A logical value indicating if mates are present in the BAM
#' file (paired-end) or not (single-end mostly).
#' @param chunkSize To safe memory the BAM file can be processes in chunks
#' of this size. Default is 3e07.
#' @param ... Further parameters that are passed to the process functions.
#' @return An RleList of coverages.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.readBAM <- function(path, indexFile = path, processFUN,
                     chromosomeList = NULL, params = NULL,
                     asMates, ...){

    ## stop if package for shifting in case of single-end data is not installed
    if(!asMates) {
        if(!requireNamespace("chipseq", quietly = TRUE)) {
            stop("The 'chipseq' package is required to estimate the fragment length of single-end data. Please install it from Bioconductor or use asMates = TRUE for paired-end data", call. = FALSE)
        }
    }

    args <- list(...)

    ## get chromosome List from BAM header and if a chromosome list is provided, filter for them
    header <- Rsamtools::scanBamHeader(path)
    chromosomeLengths <- header[[1]]$targets
    if(!is.null(chromosomeList)) {
        chromosomeLengths <- chromosomeLengths[names(chromosomeLengths) %in% chromosomeList]
    }
    
    if(is.null(params)) {
        params <- Rsamtools::ScanBamParam(what = c("pos", "qwidth"))
    }
    
    ## convert which params to GRanges
    if(length(Rsamtools::bamWhich(params)) != 0) {
        regions <- GenomicRanges::GRanges(Rsamtools::bamWhich(params))
        suppressWarnings(GenomeInfoDb::seqlengths(regions) <- chromosomeLengths[GenomeInfoDb::seqlevels(regions)])
        regions <- GenomicRanges::trim(regions)
    }
    else {
        starts <- rep(1, length(chromosomeLengths))
        ends <- chromosomeLengths
        regions <- GenomicRanges::GRanges(names(chromosomeLengths), IRanges::IRanges(starts, ends))
    }
    
    lambdaFun <- function(chromName, regions, asMates, path, indexFile, params, processFUN, args) {
        ## read data
        ##suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))
        Rsamtools::bamWhich(params) <- regions[GenomeInfoDb::seqnames(regions) == chromName,]
        if (asMates) reads <- GenomicAlignments::readGAlignmentPairs(path, index = indexFile, param = params)
        else reads <- GenomicAlignments::readGAlignments(path, index = indexFile, param = params)
        
        if(length(reads) == 0L) return(NULL)
      
        ## ensure sequence levels are correct
        GenomeInfoDb::seqlevels(reads, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(reads)
        
        ans <- do.call(processFUN, c(list(reads), args))
        return(ans)
    }
    ## run BAM reading in parallel. Check backend by bpparam().
    res <- BiocParallel::bplapply(names(chromosomeLengths), lambdaFun, regions = regions, asMates = asMates, 
                    path = path, indexFile = indexFile, params = params, 
                    processFUN = processFUN, args = args)
    ## necessary to identify null entries in list, e.g. if only a subset of 
    ## chromosomes was specified. Depending on the order of null and 
    ## non-null elements in the list the combine method will produce different
    ## results.
    idx <- !sapply(res, is.null)
    res <- do.call(c,res[idx])
    return(res)
}

#' A function to read in BED files.
#' 
#' The functions reads in BED files.
#'
#' @param path A character object indicating the path to the BED file.
#' @return A RleList of coverages.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.readBED <- function(path){}

## Process functions for chunks of data
## ====================================

#' Processing raw data and convert to integer.
#'
#' This function processes the raw data from BAM files and converts them to
#' integer type.
#' @param chunk An object of type GAlignments
#' @param paired A logical value, indicating if data is single or paired end.
#' @param chromLengths A named integer vector of chromosome lengths. 
#' @return An Rle object of integer values giving the counts for each
#' position.
#' @param center Logical value indicating if the fragments should be centered
#' or not. If centered, the fragment is reduced to only one data point,
#' which is the center of the fragment. If center = FALSE, then the complete
#' fragment is taken to compute the coverage. In single-end cases the
#' fragment size is estimated by one of the three methods: coverage,
#' correlation or SISSR. See ?chipseq::estimate.mean.fraglen.
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
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
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
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
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.getPairedCenters <- function(chunk) {
    if(length(chunk) == 0) return(GenomicRanges::GRanges())
    gr <- GenomicRanges::GRanges(chunk)

    fragmentSize <- IRanges::width(gr)
    midpoints <- (IRanges::start(gr) + IRanges::end(gr))/2
  
    IRanges::start(gr) <- IRanges::end(gr) <- midpoints

    ## filter for artefacts because of wrongly mapped reads
    maxAllowedFragSize <- 2*median(fragmentSize)
    falseFragments <- which(fragmentSize > maxAllowedFragSize)
    
    if (length(falseFragments) > 0) {
        gr <- gr[-falseFragments]
    }
    return(gr)
}

#' Extracting fragments given reads from a BAM file.
#'
#' This functions extracts the fragments from a GAligment or GAligmentPairs
#' object
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
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
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
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.getFragments <- function(chunk){
    if(length(chunk) == 0) return(GenomicRanges::GRanges())
    gr <- GenomicRanges::GRanges(chunk)

    maxAllowedFragSize <- 2*median(IRanges::width(gr))
    falseFragments <- which(IRanges::width(gr) > maxAllowedFragSize)
    if (length(falseFragments) > 0) {
        gr <- gr[-falseFragments]
    }
    return(gr)
}

## Read functions
## ==============

#' Read in raw data.
#'
#' A wrapper function to read in raw data of multiple formats (so far only
#' BAM is supported). Each list object is one chromosome.
#'
#' @param path  A character object indicating the path to the data file.
#' @param ... Further arguments to be passed to the low-level read functions.
#' @return An RleList of coverages.
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
.readRawData <- function(path, ...) {
    elements <- strsplit(path,"\\.")[[1]]
    suffix <- elements[length(elements)]
    res <- NULL
    if (suffix == "bam") {
        res <- .readBAM(path = path, ...)
    }
    if (suffix == "bed") {
        res <- .readBED()
    }
    return(res)
}
    
#' Read and process raw data
#'
#' This is a convenience function for reading and parsing raw data from
#' multiple formats (at the moment only BAM is supported) and process to
#' GenomicTiles object for further analysis.
#'
#' @param path A character object indicating the path to the source file
#' @param chromosomeList A character vector of chromosomes.
#' @param asMates Is the data paired end?
#' @param processFUN A function specifying how to process raw data. There
#' are three restriction. First, the first argument of the function must take
#' a GAligments object, as this is the output of readGAligments(Pairs).
#' Second, the output must be an Rle vector of sorted positions, e.g.
#' 345, 345, 347, 349, ... Third, each processed Rle vector, must contain
#' a metadata list with the elements 'start', 'end' and 'chrom', specifying
#' at which position the vector starts and ends and on which chromosome it
#' lies. 
#' @param ... Further parameters that can be passed to low-level functions.
#' Mostly to pass arguments to custom process functions. In case the default
#' process functions are used the most interesting parameters might be
#' from ?chipseq::estimate.mean.fraglen for single-end data.
#' @return An RleList of coverages per chromosome.
#' @examples
#' \dontrun{
#' path <- "path/to/file.bam"
#' gt <- readData(path)
#' header <- Rsamtools::scanBamHeader(pathSingle)
#' chromosomeList <- header[[1]]$targets
#' gt <- readData(path,chromosomeList)
#' }
#' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
.readData <- function(path, chromosomeList = NULL, asMates = TRUE, processFUN = NULL, ...) {

    header <- Rsamtools::scanBamHeader(path)
    chromosomeLengths <- header[[1]]$targets
    
    if(!is.null(chromosomeList)) {
        chromosomeLengths <- chromosomeLengths[names(chromosomeLengths) %in% chromosomeList]
    }

    data <- do.call(.readRawData,
                    c(list(path = path, chromosomeList = chromosomeList,
                           processFUN = processFUN, asMates = asMates),
                      list(...)))
    seqlevels(data) <- seqlevelsInUse(data)
    res <- IRanges::coverage(data)
    #res <- .combineTrack(data, chromosomeLengths)
    
    return(res)
}

## Mappability (deprecated for now)
## ===============================

## #' Read mappability data.
## #'
## #' The function reads in a BAM file and returns a list of logical vectors for each chromosome,
## #' indicating if the position is mappable or not. This function is not intended for direct use.
## #' However one might parse this data differently, than the package does. Therefore the
## #' function is made available.
## #'
## #' @param path A character object indicating the path to the BAM file.
## #' @param ... Any other parameters supplying to \code{\link{readRawData}}.
## #' @return An RleList consisting of logical-Rle objects, where each objects has the length
## #' length(chromosome) representing if a position on the given chromosome is mappable or not.
## #' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
## readMapBAM <- function(path, chromosomeList = NULL, ...) {

##     if(is.null(chromosomeList)) {
##         header <- Rsamtools::scanBamHeader(path)
##         chromosomeList <- header[[1]]$targets
##     }
    
##     temp <- readRawData(path, chromosomeList = chromosomeList, processFUN = .processMappabilityChunks, ...)
##     res <- .combineTrack(temp, chromosomeList, init = FALSE)
##     return(res)
## }

## #' Processing raw data and convert to logical.
## #'
## #' This function processes the data read from BAM files and converts them to logical type.
## #' @param chunk An object of type GAlignments
## #' @return An Rle object of logical values indicating for each position, if it is mappabale or not.
## #' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
## .processMappabilityChunks <- function(chunk) {
##     centeredFragments <- .centerFragments(chunk, asMates = TRUE)
##     start <- min(centeredFragments)
##     end <- max(centeredFragments)
##     chrom <- as.character(unique(seqnames(chunk)))
##     track <- Rle(FALSE, end - start + 1)
##     track[(runValue(centeredFragments) - start + 1)] <- TRUE
##     track@metadata <- list(start = start, end = end, chromosome = chrom)
##     return(track)
## }


## #' A function to parse centered BAM files to a GRanges object of mappable regions.
## #'
## #' The function takes in a logical-RleList indicating the mappability of each position in the genome and returns
## #' a GRanges object of mappable regions. This function is not intended to be used directly.
## #'
## #' @param obj A logical RleList, where each list object represents one chromosome
## #' @param winsize The size for the moving window. Should be odd, one will be added if this isn't the case.
## #' Default: 101
## #' @param mincov The minimal coverage of TRUE values within the moving window as a fraction. Every position
## #' above this value will be regarded as mappable.
## #' @param minlength The minimal length of consecutive mappabale positions that form one mappable region.
## #' Every region below that length will be marked unmappable.
## #' @return A GRanges object of mappable regions
## #' @details A moving window of size 'winsize' is used to compute the fraction of mappable position in the
## #' genome. Each position above the 'mincov' value is regarded as mappable. Consecutive mappable positions
## #' are gathered to form a region. if the region is longer than 'minlength' it forms a mappable region.
## #' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
## parseMappability <- function(obj, winsize = 101, mincov = 0.5, minlength = 1000) {
##     if (winsize %% 2 == 0) winsize <- winsize+1

##     chromosomeList <- names(obj)
##     res <- bplapply(chromosomeList, function(chromosome) {
##         require(genoGAM, quietly = TRUE)
##         binary_mappability <- obj[[chromosome]]*1
##         mappability_MA <- runmean(binary_mappability, k = winsize, endrule = "constant")
##         boolean_mappability <- Rle(ifelse(mappability_MA < mincov, FALSE, TRUE))
        
##         ## Set TRUE stretches of < MIN_LENGTH to FALSE:     
##         indx <- runValue(boolean_mappability) == TRUE & runLength(boolean_mappability) < minlength
##         runValue(boolean_mappability)[indx] <- FALSE
        
##         mappability_region <- GRanges(seqnames = chromosome,
##                                       ranges = IRanges(start = start(boolean_mappability),
##                                           end = end(boolean_mappability)),
##                                       mapable = runValue(boolean_mappability))
        
##         return(mappability_region)
##     })
##     names(res) <- chromosomeList
##     return(GRangesList(res))
## }

## #' Read and convert raw mappability data
## #'
## #' This is a convenience function for reading and parsing raw data from multiple formats (at the moment
## #' only BAM is supported) and converting to a list containing the mappability track and mappable region
## #'
## #' @param path A character object indicating the path to the source file
## #' @param winsize The size for the moving window. Should be odd, one will be added if this isn't the case.
## #' Default: 101
## #' @param mincov The minimal coverage of TRUE values within the moving window as a fraction. Every position
## #' above this value will be regarded as mappable.
## #' @param minlength The minimal length of consecutive mappabale positions that form one mappable region.
## #' Every region below that length will be marked unmappable.
## #' @param ... Further parameters that can be passed to \code{\link{readMapBAM}}.
## #' @return A GRanges object of mappable regions
## #' @examples
## #' \dontrun{
## #' path <- "path/to/mappability.bam"
## #' map <- readMappability(path)
## #' map <- readData(path,asMates=FALSE,readByChromosome = TRUE)
## #' mappabilityTrack <- map$track
## #' mappabilityRegions <- map$regions
## #' }
## #' @author Georg Stricker \email{stricker@@genzentrum.lmu.de}
## readMappability <- function(path, winsize = 101, mincov = 0.5, minlength = 1000, ...) {

##     data <- readMapBAM(path, ...)
##     mapRegions <- parseMappability(data, winsize = winsize, mincov = mincov, minlength = minlength)
##     res <- list(track = data, regions = mapRegions)
    
##     return(res)
## }

################## Tests
# ======================

library(SummarizedExperiment)
library(HDF5Array)
FOLDER <- "/s/project/coreProm/Michi/thornton/align_STAR"
CONFIG <- "~/workspace/analysis/diffBinding/config.txt"
config <- data.table::fread(CONFIG)
rhdf5::H5close()

BiocParallel::register(BiocParallel::SnowParam(workers=4))
## first approach: load all data in memory and then write to HDF5
## might be memory consuming, especially for human data
# ======================

fullRead <- function(config, folder) {
    header <- Rsamtools::scanBamHeader(file.path(folder, config$file[1]))
    chromosomeLengths <- header[[1]]$targets
    len <- sum(chromosomeLengths)
    sampleList <- vector("list", nrow(config))

    for(ii in 1:nrow(config)) {
        path <- file.path(folder, config$file[ii])
        reads <- GenomicAlignments::readGAlignments(path)
        counts <- coverage(reads)
        counts <- unlist(counts)
        zeros <- len - length(counts)
        sampleList[[ii]] <- c(counts, rep(0, zeros))
    }
    names(sampleList) <- config$ID
    h5 <- HDF5Array(DelayedArray(DataFrame(sampleList)))
    
    se <- SummarizedExperiment(h5)
    rowRanges(se) <- GPos(GRanges(names(chromosomeLengths), IRanges(rep(1, length(chromosomeLengths)), chromosomeLengths)))
    return(se)
}

## second approach: load only one sample to memory, wrire to HDF5 and then cbind all other samples
# =========================

partialRead <- function(config, folder) {
    h5 <- NULL
    header <- Rsamtools::scanBamHeader(file.path(folder, config$file[1]))
    chromosomeLengths <- header[[1]]$targets
    len <- sum(as.numeric(chromosomeLengths))
    
    for(ii in 1:nrow(config)) {
        path <- file.path(folder, config$file[ii])
        reads <- GenomicAlignments::readGAlignments(path)
        counts <- coverage(reads)
        counts <- unlist(counts)
        zeros <- len - length(counts)
        counts <- c(counts, rep(0, zeros))
        df <- DataFrame(sample1 = counts)
        names(df) <- config$ID[ii]
        da <- DelayedArray(df)
        if(ii == 1) {
            h5 <- HDF5Array(da)
        }
        else {
            temp <- HDF5Array(da)
            h5 <- cbind(h5, temp)
        }
    }
    se <- SummarizedExperiment(h5)
    rowRanges(se) <- GPos(GRanges(names(chromosomeLengths), IRanges(rep(1, length(chromosomeLengths)), chromosomeLengths)))
    return(se)
}

## Comparison
library(profvis)

profvis({
    df <- partialRead(config, FOLDER)
    df <- fullRead(config, FOLDER)
})

## no real difference on yeast. A little less memory usage for partial, but a little bit slower than full.

## Comparison for human
# =====================
## For human the comparison is between partial (chromosome wise) and full (whole sample).
## Reading first the entire data frame and writing then is not memory efficient
## One entire genome of human is around 24GB = not efficient to read and write it entirely
FOLDER <- "/s/project/coreProm/data/hg38"
CONFIG <- "~/workspace/analysis/peakcalls/config_cebpb.txt"
config <- data.table::fread(CONFIG)
rhdf5::H5close()

df <- partialRead(config, FOLDER)
