###############################
## p-value computation
###############################

#' computeSignificance
#'
#' The function computes positionwise p-values
#' 
#' @param gg A GenoGAM object.
#' @param log.p Should values be returned in log scale?
#' @return An updated GenoGAM object, where the pvalue slot is added.
#' @details Note, that in case the data is stored in HDF5 format, the pvalue 'group' is added on hard drive.
#' That is, unlike any other function in R, where the input object is not changed, it actually
#' is in this case. If one wishes to have HDF5 data without the pvalue 'group', one has to
#' backup the HDF5 files prior to computation or delete them after with \code{rhdf5::h5delete}
#' @examples
#' ## make test GenoGAM
#' gg <- makeTestGenoGAM()
#' ## compute pvalues
#' computeSignificance(gg)
#' pvalue(gg)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
computeSignificance <- function(gg, log.p = FALSE) {

    if(log.p) {
        futile.logger::flog.info("Computing positionwise log p-values")
    }
    else {
        futile.logger::flog.info("Computing positionwise p-values")
    }

    res <- .pvals(gg, log.p)
    
    futile.logger::flog.info("DONE")
    return(res)
}

.pvals <- function(gg, log.p = FALSE) {
    hdf5 <- is.HDF5(gg)
    split <- is(gg, "GenoGAMList")

    if(split){
        
        if(hdf5) {
            futile.logger::flog.debug("Data is in HDF5 format and split by chromosome")
            res <- .pvals_hdf5_split(gg, log.p)
        }
        else {
            futile.logger::flog.debug("Data is split by chromosome")
            res <- .pvals_split(gg, log.p)
        }
    }
    else {
        if(hdf5) {
            futile.logger::flog.debug("Data is in HDF5 format")
            res <- .pvals_hdf5(gg, log.p)
        }
        else {
            res <- .pvals_default(gg, log.p)
        }
    }
    return(res)
}

## compute pvalue for GenoGAM object without HDF5 and no split
.pvals_default <- function(gg, log.p = FALSE) {
    tracks <- colnames(gg)
    pvals <- lapply(tracks, function(id) {
        ## base <- median(fits(gg)[[id]]) ## get the logarithmic base level of the background
        base <- 0
        ## futile.logger::flog.debug(paste("The base level for track", id, "computed as:", base))
        
        res <- 2*pnorm(base, mean = abs(fits(gg)[[id]]), sd = se(gg)[[id]], log.p = log.p)
        return(res)
    })
    df <- S4Vectors::DataFrame(data.frame(pvals))
    names(df) <- tracks

    ## add pvals to assay fields. SummarizedExperiment has some control functions
    ## guarding the structure of the data, so we have to directly plug it in to
    ## avoid copying
    gg@assays@data@listData$pvalue <- df
    return(gg)
}

## compute pvalue for GenoGAM object with HDF5 and no split
.pvals_hdf5 <- function(gg, log.p = FALSE) {
    tracks <- colnames(gg)

    ## make HDF5 group in existing HDF5 file for p-values
    dims <- dim(SummarizedExperiment::assay(gg))
    
    seedFile <- .get_seed(SummarizedExperiment::assay(gg))
    h5name <- "pvalue"
    
    ## check if h5name already exists
    h5content <- rhdf5::h5ls(seedFile)
    if(h5name %in% h5content$name) {
        futile.logger::flog.warn(paste0("HDF5 matrix '", h5name, "' already exists. Overwriting"))
        rhdf5::h5delete(seedFile, h5name)
    }
    
    h5file <- rhdf5::H5Fopen(seedFile)
    chunks <- HDF5Array::getHDF5DumpChunkDim(dims)
   
    .createH5Dataset(h5file, name = h5name, settings = gg@settings,
                     d = dims, chunk = chunks, type = "H5T_IEEE_F32LE")
    
    ## close file to save initial zero matrix
    rhdf5::H5Fclose(h5file)
    
    ## reopen to write again
    h5file <- rhdf5::H5Fopen(seedFile)
    pvals <- lapply(1:length(tracks), function(id) {
        sp_fits <- fits(gg)[,id, drop = FALSE]
        sp_ses <- se(gg)[,id, drop = FALSE]

        base <- integer(length(tracks))
        futile.logger::flog.debug(paste("The base level for track", tracks[id], "computed as:", base))
        
        .hdf5_block_pval(x = sp_fits, se = sp_ses, base = base,
                         h5file = h5file, name = h5name, id = id,
                         log.p = log.p)
    })

    ## add pvals to assay fields. SummarizedExperiment has some control functions
    ## guarding the structure of the data, so we have to directly plug it in to
    ## avoid copying
    h5 <- HDF5Array::HDF5Array(seedFile, h5name)
    gg@assays@data@listData$pvalue <- h5

    ## close again
    rhdf5::H5Fclose(h5file)

    return(gg)
}

## compute pvalue for GenoGAM object with HDF5 and split
.pvals_hdf5_split <- function(gg, log.p = FALSE) {
    tracks <- colnames(gg)
    sp_fits <- fits(gg)
    sp_ses <- se(gg)
    chroms <- names(sp_fits)

    ## set base level
    base <- integer(length(tracks))
    if(futile.logger::flog.threshold() == "DEBUG") {
        for(ii in 1:length(base)) {
            futile.logger::flog.debug(paste("The base level for track", tracks[ii], "computed as:", base[[ii]]))
        }
    }

    pvals_list <- lapply(chroms, function(y) {
        futile.logger::flog.debug(paste("Computing p-values for chromosome", y))
        ## make HDF5 group in existing HDF5 file for p-values
        dims <- dim(SummarizedExperiment::assay(gg@data[[y]]))
    
        seedFile <- .get_seed(SummarizedExperiment::assay(gg@data[[y]]))
        h5name <- "pvalue"

        ## check if h5name already exists
        h5content <- rhdf5::h5ls(seedFile)
        if(h5name %in% h5content$name) {
            futile.logger::flog.warn(paste0("HDF5 matrix '", h5name, "' in chromosome ", y, " already exists. Overwriting"))
            rhdf5::h5delete(seedFile, h5name)
        }
        
        h5file <- rhdf5::H5Fopen(seedFile)
        chunks <- .get_chunkdim(SummarizedExperiment::assay(gg@data[[y]]))
        
        .createH5Dataset(h5file, name = h5name, settings = gg@settings, d = dims,
                         chunk = chunks, type = "H5T_IEEE_F32LE")
        ## close file to save initial zero matrix
        rhdf5::H5Fclose(h5file)
        
        ## reopen to write again
        h5file <- rhdf5::H5Fopen(seedFile)
        pvals <- lapply(1:length(tracks), function(id) { 
            .hdf5_block_pval(x = sp_fits[[y]][,id, drop = FALSE],
                             se = sp_ses[[y]][,id, drop = FALSE],
                             base = base[[id]], h5file = h5file,
                             name = h5name, id = id, log.p = log.p)
        })

        ## add pvals to assay fields. SummarizedExperiment has some control functions
        ## guarding the structure of the data, so we have to directly plug it in to
        ## avoid copying
        h5 <- HDF5Array::HDF5Array(seedFile, h5name)

        ## close again
        rhdf5::H5Fclose(h5file)
	return(h5)
    })

    for(ii in 1:length(chrom)) {
        gg@data[[chrom[ii]]]@assays@data@listData$pvalue <- pvals_list[[ii]]
    }

    return(gg)
}

.hdf5_block_pval <- function(x, se, base, h5file, name, id, log.p = FALSE) {

    grid <- .makeGrid(x)
    nblock <- length(grid)

    ## now compute pvalues
    for (b in seq_len(nblock)) {
        subgrid <- grid[[b]]
        xblock <- DelayedArray:::read_block(x, subgrid)
        seblock <- DelayedArray:::read_block(se, subgrid)
        xtmp <- as.vector(xblock, mode = "numeric")
        setmp <- as.vector(seblock, mode = "numeric")

        res <- 2*pnorm(base, mean = abs(xtmp), sd = setmp, log.p = log.p)
        idx <- IRanges::ranges(subgrid)[1]
        rhdf5::h5write(as.matrix(res, ncol = 1), file = h5file, name = name,
                       index = list(start(idx):end(idx), id))
    }
    invisible(NULL)
}

## compute pvalue for GenoGAM object without HDF5, but with Split
.pvals_split <- function(gg, log.p = FALSE) {
    tracks <- colnames(gg)
    chroms <- names(SummarizedExperiment::assays(gg))
    len <- length(chroms)

    ## set base
    base <- list(0, 0)
    names(base) <- tracks

    ## now compute p-values
    pvals_list <- lapply(chroms, function(y) {
        futile.logger::flog.debug(paste("Computing p-values for chromosome", y))
        pvals <- lapply(tracks, function(id) {
            res <- 2*pnorm(base[[id]], mean = abs(fits(gg)[[y]][[id]]), sd = se(gg)[[y]][[id]], log.p = log.p)
            return(res)
        })
        df <- DataFrame(data.frame(pvals))
        names(df) <- tracks

        ## add pvals to assay fields. SummarizedExperiment has some control functions
        ## guarding the structure of the data, so we have to directly plug it in to
        ## avoid copying
        return(df)
    })
    for(ii in 1:length(chroms)) {
        gg@data[[chroms[ii]]]@assays@data@listData$pvalue <- pvals_list[[ii]]
    }
    return(gg)
}
