#' make HDF5 name for the different files to be written
#' @noRd
.makeFilename <- function(dir, name, seed = NULL, split = " ", id = FALSE) {
    if(!is.null(seed)) {
        ident <- strsplit(seed, split = split)[[1]]
        suffix <- ident[length(ident)]
        f <- paste0(name, "_", suffix)
    }
    else {
        if(id) {
            date <- strsplit(as.character(Sys.time()), " ")[[1]][1]
            suffix <- paste0(gsub("-", "", date), ".h5")
            f <- paste0(name, "_", suffix)
        }
        else {
            f <- name
        }
    }

    h5 <- file.path(dir, f)
    return(h5)
}

#' Creates an HDF5 dataset on hard drive from name and with specified settings
#' @noRd
.createH5Dataset <- function(h5file, name, settings, d, chunk, type) {
    
    level <- slot(settings, "hdf5Control")$level
    
    if(!rhdf5::h5createDataset(h5file, name, d, d,
                           H5type = type, chunk = chunk,
                           level = level)) {
        stop("Couldn't create HDF5 dataset.")
    }
    
    invisible(NULL)
}

##' Function to write DataFrame of pre-processed counts to HDF5
##' @noRd
.writeToHDF5 <- function(df, file, chunks, name = file, settings, simple = FALSE, type = "H5T_STD_I32LE") {

    ## create dir if not present
    dir <- slot(settings, "hdf5Control")$dir
    if(!dir.exists(dir)) {
        futile.logger::flog.info(paste("HDF5 directory created at:", dir))
        dir.create(dir)
    }

    futile.logger::flog.info(paste("Writing", name, "to HDF5"))

    ## make HDF5 chunks
    d <- dim(df)
    h5chunk <- c(max(width(chunks)), ncol(df))

    ## create HDF5 file
    h5file <- .createH5File(dir, file, id = !simple)

    ## initialize HDF5 dataset
    .createH5Dataset(h5file$pointer, name, settings, d,
                     type = type, chunk = h5chunk)
    
    ## split to writeable chunks
    rle <- extractList(df, ranges(chunks))

    ## write to file
    sapply(1:length(rle), function(y) {
        rhdf5::h5write(as.matrix(as.data.frame(rle[[y]])), file = h5file$pointer, name = name,
                       index = list(start(chunks[y,]):end(chunks[y,]), 1:d[2]))
    })
    rhdf5::H5Fclose(h5file$pointer)
    
    ## link to DelayedArray
    h5 <- HDF5Array::HDF5Array(h5file$file, name = name)
    
    futile.logger::flog.info(paste(name, "written"))
    
    return(h5)
}


#' create HDF5 file
#' @noRd
.createH5File <- function(dir, name, seed = NULL, split = " ", id = FALSE) {

    f <- .makeFilename(dir, name, seed = seed, split = split, id = id)
    if(file.exists(f)) {
        warning(paste("File:", f, "exists. Overwriting"))
        rhdf5::H5close()
        file.remove(f)
    }
    h5file <- rhdf5::H5Fcreate(f)
    
    return(list(pointer = h5file, file = f))
}

#' initialize HDF5 dataset with specified dimensions and chunk
#' @noRd
.createH5DF <- function(seed, settings, d, chunk = NULL, what = c("fits", "coefs", "sumMatrix")) {
    what <- match.arg(what)

    dir <- slot(settings, "hdf5Control")$dir
    
    if(what == "fits") {
        h5file <- .createH5File(dir, name = what, seed = seed, split = "/")
        .createH5Dataset(h5file$pointer, name = "fits", settings = settings, d = d,
                         chunk = chunk, type = "H5T_IEEE_F32LE")
        .createH5Dataset(h5file$pointer, name = "ses", settings = settings, d = d,
                         chunk = chunk, type = "H5T_IEEE_F32LE")
    }
    if(what == "coefs") {
        h5file <- .createH5File(dir, name = what, seed = seed, split = "_")
        .createH5Dataset(h5file$pointer, name = "coefs", settings = settings, d = d,
                         chunk = chunk, type = "H5T_IEEE_F32LE")
    }
    if(what == "sumMatrix") {
        h5file <- .createH5File(dir, name = seed, id = TRUE)
        .createH5Dataset(h5file$pointer, name = seed, settings = settings, d = d,
                         chunk = chunk, type = "H5T_STD_I32LE")
    }
    rhdf5::H5Fclose(h5file$pointer)

    return(h5file$file)
}

## Queue mechanism
## -------------------------------------------------------------------------

#' initialized queue
#' @noRd
.init_Queue <- function(file){
    if(is.null(file)) {
        return(NULL)
    }
    
    dir <- strsplit(file, split = "\\.")[[1]][1]
    if(dir.exists(dir)) {
        warning(paste("Directory:", dir, "exists. Overwriting"))
        unlink(dir, recursive = TRUE)
    }
    dir.create(dir)
    return(dir)
}

#' push to queue
#' @noRd
.queue <- function(dir) {
    qid <- as.integer(Sys.time())
    pid <- Sys.getpid()
    f <- paste(qid, pid, sep = "_")
    if(!file.create(file.path(dir,f))) {
        stop(paste("Could not queue. Check if directory folder:", dir, "is correct or has writing permissions"))
    }
    Sys.sleep(0.5)
    return(f)
}

#' pop from queue
#' @noRd
.unqueue <- function(qid, dir) {
    f <- file.path(dir, qid)
    if(!file.remove(f)) {
        stop(paste("Could not unqueue. Check if file:", f, "exists the writing permissions"))
    }
    invisible()
}

#' close queue
#' @noRd
.end_Queue <- function(dir) {
    if(!is.null(dir)) {
        unlink(dir, recursive = TRUE)
    }
    
    invisible()
}

#' Helper function to make a grid for HDF5 blocks.
#' Grid is used to move from block to block and apply functions
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.makeGrid <- function(x) {
    max_block_size <- DelayedArray:::getAutoBlockLength(type(x))
    xsize <- nrow(x)
    mult <- as.integer(ceiling(xsize/max_block_size))
    optimal_block_size <- as.integer(xsize/mult)
    grid <- DelayedArray::blockGrid(x, optimal_block_size)
    
    return(grid)
}

#' The essential block_APPLY function that controls functions to be applied
#' on a block which work in a map-reduce manner
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.block_APPLY <- function(x, fun, ...) {

    grid <- .makeGrid(x)
    nblock <- length(grid)

    qn <- numeric(nblock)
    for (b in seq_len(nblock)) {
        ## take block and apply function
        block <- DelayedArray:::read_block(x, grid[[b]])
        qn[b] <- fun$APPLY(block, ...)
    }
    res <- fun$COMBINE(qn)
    
    return(res)
}

#' A helper function to get the filepath of an HDF5 file
#' @noRd
.get_seed <- function(x) {
    if("filepath" %in% slotNames(x)) {
        return(x@filepath)
    }
    else {
        .get_seed(x@seed)
    }
}

#' A helper function to get the filepath of an HDF5 file
#' @noRd
.get_chunkdim <- function(x) {
    if("chunkdim" %in% slotNames(x)) {
        return(x@chunkdim)
    }
    else {
        .get_chunkdim(x@seed)
    }
}

### The folowing is modified from HDF5Array package. Thanks to Herve Pagès.
## NOT USED as locking does not work for too many parallel workers

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A simple mechanism to lock/unlock a file so processes can get temporary
### exclusive access to it
###

## ##' creates the lock file
## .locked_path <- function(filepath)
## {
##     if (!isSingleString(filepath) || filepath == "") {
##         stop("'filepath' must be a single non-empty string")
##     }
##     dir <- dirname(filepath)
##     return(file.path(dir, ".lock"))
## }

## ##' locks the file by creating an empty ".lock" file
## .lock_file <- function(filepath)
## {
##     locked_path <- .locked_path(filepath)
    
##     ## Must wait if the file is already locked.
##     while (file.exists(locked_path)) {
##         ## Sys.sleep(runif(1, 0.1, 1))
##         Sys.sleep(2)
##     }

##     if(!file.create(locked_path)) {
##         stop("failed to lock '", filepath, "' file")
##     }
##     return(locked_path)
## }

## ##' unlocks file by removing the empty ".lock" file
## .unlock_file <- function(lock)
## {
##     if (!file.remove(lock)) {
##         stop("failed to unlock file")
##     }
## }
