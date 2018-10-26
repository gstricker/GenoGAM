#################################
## The main modelling function ##
#################################

#' genogam
#'
#' The main modelling function.
#' @param ggd The GenoGAMDataSet object to be fitted
#' @param lambda The penalization parameter. If NULL (default) estimated by
#' cross validation. So far only one parameter for all splines is supported.
#' @param theta The global overdispersion parameter. If NULL (default) estimated
#' by cross validation.
#' @param family The distribution to be used. So far only Negative-Binomial (nb)
#' is supported.
#' @param eps The factor for additional first-order regularization. This should be
#' zero (default) in most cases. It can be useful for sparse data with many
#' zero-counts regions or very low coverage. In this cases it is advised to use
#' a small factor like 0.01, which would penalize those regions but not the ones
#' with higher coverage. See Wood S., Generalized Additive Models (2006) for more.
#' @param kfolds The number of folds for cross validation
#' @param intervalSize The size of the hold-out intervals in cross validation.
#' If replicates are present, this can easily be increased to twice the fragment
#' size to capture more of the local correlation. If no replicates are present,
#' keep the number low to avoid heavy interpolation (default).
#' @param regions How many regions should be used in cross validation? The number
#' is an upper limit. It is usually corrected down, such that the total number of
#' models computed during cross validation does not exceed the total number of
#' models to compute for the entire genome. This is usually the case for small
#' organisms such as yeast.
#' @param order The order of the B-spline basis, which is order + 2, where 0
#' is the lowest order. Thus order = 2 is equivalent to cubic order (= 3).
#' @param m The order of penalization. Thus m = 2 penalizes the second differences.
#' @return The fit as a GenoGAM object
#' @examples
#' ggd <- makeTestGenoGAMDataSet(sim = TRUE)
#' res <- genogam(ggd, lambda = 266.8368, theta = 2.415738)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname fitGenoGAM
#' @export
genogam <- function(ggd, lambda = NULL, theta = NULL, family = "nb", eps = 0,
                    kfolds = 10, intervalSize = 20, regions = 20, order = 2,
                    m = 2) {

    futile.logger::flog.info("Initializing the model")

    input <- paste0("Initializing the model with the following parameters:\n",
                    "  Lambda: ", lambda, "\n",
                    "  Theta: ", theta, "\n",
                    "  Family: ", family, "\n",
                    "  Epsilon: ", eps, "\n",
                    "  Number of folds: ", kfolds, "\n",
                    "  Interval size: ", intervalSize, "\n",
                    "  Number of regions: ", regions, "\n",
                    "  B-spline order: ", order, "\n",
                    "  Penalization order: ", m, "\n",
                    "  GenoGAMDataSet object: ")
    futile.logger::flog.debug(input)
    futile.logger::flog.debug(show(ggd))
    futile.logger::flog.debug(show(slot(ggd, "settings")))
    
    settings <- slot(ggd, "settings")
    check <- .checkSettings(ggd)
    coords <- .getCoordinates(ggd)
    chunks <- .getChunkCoords(coords)

    futile.logger::flog.info("Done")
    ## Cross Validation
    cv <- FALSE
    regionSize <- slot(settings, "dataControl")$regionSize
    bpknots <- slot(settings, "dataControl")$bpknots

    ggs <- setupGenoGAM(ggd, lambda = lambda, theta = theta, family = family, 
                        eps = eps, bpknots = bpknots, order = order,
                        penorder = m, control = slot(settings, "estimControl"))
    
    if(is.null(lambda) | is.null(theta)) {
        futile.logger::flog.info("Estimating hyperparameters")

        ## check correct region size
        if(regionSize < 2*getOverhangSize(ggd)) {
            futile.logger::flog.warn("region size too small for cross validation, set to twice the overhang size.")
            regionSize <- 2*getOverhangSize(ggd)
        }

        ## backup original tile index
        index_backup <- getIndex(ggd)

        ## make new tile index
        metadata(slot(ggd, "index"))$chunkSize <- regionSize
        newTiles <- .makeTiles(tileSettings(ggd))
        slot(ggd, "index") <- newTiles
        new_ggs <- setupGenoGAM(ggd, lambda = lambda, theta = theta, family = family, 
                                eps = eps, bpknots = bpknots, order = order,
                                penorder = m, control = slot(settings, "estimControl"))
        new_coords <- .getCoordinates(ggd)
        
        ## get the tile ids for CV
        sumMatrix <- getCountMatrix(ggd)
        
        ncv <- regions
        
        if(ncv < length(new_coords)) {
            if(sum(sapply(colData(ggd), sum)) == nrow(colData(ggd)) |
               nrow(sumMatrix) < regions) {
                if(is(sumMatrix, "DelayedMatrix")) {
                    rsums <- DelayedArray::rowSums(sumMatrix)
                }
                else {
                    rsums <- rowSums(sumMatrix)
                }
                ids <- order(rsums, decreasing = TRUE)[1:ncv]
            }
            else {
                futile.logger::flog.debug("Performing DESeq2 analysis to find regions with highest fold change")
                pvals <- suppressMessages(suppressWarnings(.deseq(sumMatrix, colData(ggd))))
                ids <- order(pvals)[1:ncv]
            }
        }
        else ids <- 1:length(new_coords)

        control <- slot(settings, "optimControl")
        futile.logger::flog.debug(paste("Selected the following regions:", paste(ids, collapse = ",")))
        
        params <- .doCrossValidation(ggd, setup = new_ggs, coords = new_coords, 
                                     id = ids, folds = kfolds, 
                                     intervalSize = intervalSize,
                                     fn = .loglik, 
                                     method = slot(settings, "optimMethod"),
                                     control = control) 

        names(params) <- c("lambda", "theta")
        slot(ggs, "params")$lambda <- params["lambda"]
        slot(ggs, "params")$theta <- params["theta"]

        ## set back the original tile index and delete other objects
        slot(ggd, "index") <- index_backup
        rm("new_ggs", "new_coords")
        gc()
        
        cv <- TRUE
        futile.logger::flog.info("Done")
    }

    ## specify ids
    ids <- 1:length(coords)
    ## make chunk coordinates relative
    s <- start(chunks) %% start(coords) + 1
    e <- s + width(chunks) - 1
    relativeChunks <- IRanges::IRanges(s, e)

    ## prepare some slots
    modelParams <- slot(ggs, "params")
    modelParams$cv <- cv
    modelParams$bpknots <- bpknots
    modelParams$order <- order
    modelParams$penorder <- m
    if(cv) {
        modelParams$kfolds <- kfolds
        modelParams$intervalSize <- intervalSize
        modelParams$regions <- ncv
    }
    modelParams <- c(modelParams, tileSettings(ggd))
    modelParams$check <- NULL

    futile.logger::flog.info("Fitting model")
    
    ## ================
    ## start model
    ## ================

    ## For a dataset split by chromosomes
    ## ----------------------------------
    
    if(is(ggd, "GenoGAMDataSetList")) {
        identifier <- slot(ggd, "id")$id
        indx <- getIndex(ggd)

        ## And with HDF5 backend
        ## ---------------------
        if(slot(ggd, "hdf5")) {
            ## create HDF5 file for coefficients
            ## the dimension of the matrix for coefsfile (betas * tile width)
            nfun <- length(.getVars(design(ggd), type = "covar"))
            nbetas <- dim(slot(ggs, "designMatrix"))[2]
            d <- c(nbetas * nfun, length(getIndex(ggd)))

            ## create Coefs file
            seedFile <- .get_seed(assay(ggd)[[1]])
            chunk <- c(nbetas * nfun, 1)
            coefsFile <- .createH5DF(seedFile, settings, d, chunk, what = "coefs")
        }

        ## lapply by identifier, i.e. chromosome
        selist <- lapply(identifier, function(y) {
            ## get correct ids
            subids  <- ids[as.vector(GenomeInfoDb::seqnames(indx) %in% y)]
            rr <- rowRanges(ggd)[[y]]
            
            ## compute model for one identifier, i.e chromosome
            futile.logger::flog.info(paste("Fitting", y))

            ## And with HDF5 backend
            ## ---------------------
            
            if(slot(ggd, "hdf5")) {               
                ## the dimension of the matrix for given chromosome
                d <- c(length(rr), nfun)
                ## create datasets
                seedFile <- .get_seed(assay(ggd)[[y]])
                h5file <- .createH5DF(seedFile, settings, d, chunk = c(getChunkSize(ggd), nfun))

                ## create chunks coordinates for given chromosome
                chromIndex <- getIndex(ggd)[seqnames(getIndex(ggd)) == y,]
                seqlevels(chromIndex, pruning.mode = "coarse") <- seqlevelsInUse(chromIndex)
                chunks <- .getChunkCoords(chromIndex)
            }

            ## No HDF5 backend
            ## ---------------------
            
            else {
                h5file <- NULL
                coefsFile <- NULL
            }

            ## initialize queue and start fitting in parallel
            qdir <- .init_Queue(h5file)
            res <- BiocParallel::bplapply(subids, .fitGenoGAM, 
                                          data = ggd, init = ggs, coords = coords,
                                          relativeChunks = relativeChunks, h5file = h5file,
                                          chunks = chunks, coefsFile = coefsFile,
                                          qdir = qdir)
            ## remove temporary queue folder and close files
            .end_Queue(qdir)
            
            futile.logger::flog.info(paste(y, "Done"))
      
            ## assemble results
            futile.logger::flog.info("Processing Fits")

            ## With HDF5 backend
            ## ---------------------
            
            if(slot(ggd, "hdf5")) {

                ## make names, as HDF5 does not store them
                colDataNames <- .makeNames(design(ggd))
                df <- DataFrame(matrix(,length(colDataNames),0))
                rownames(df) <- colDataNames

                ## TODO: In Hdf5 write function need to attach log entry to HDF5 dumpLog
                h5fits <- HDF5Array::HDF5Array(h5file, name = "fits")
                h5ses <- HDF5Array::HDF5Array(h5file, name = "ses")
                se <- SummarizedExperiment::SummarizedExperiment(rowRanges = rr,
                                                                 assays = list(fits = h5fits,
                                                                               se = h5ses))
                colData(se) <- df
                return(se)
            }

            ## No HDF5 backend
            ## ---------------------
            
            else {
                combinedFits <- .transformResults(res, relativeChunks, what = "fits")
                combinedSEs <- .transformResults(res, relativeChunks, what = "se")
                se <- SummarizedExperiment::SummarizedExperiment(rowRanges = rr,
                                                                 assays = list(fits = combinedFits,
                                                                               se = combinedSEs))
                ## extract coefficients from fits
                coefs <- sapply(res, function(y) {
                    slot(y, "beta")
                })
                return(list(se = se, coefs = coefs))
            }
        })

        if(slot(ggd, "hdf5")) {
            names(selist) <- names(assay(ggd))
            coefs <- HDF5Array::HDF5Array(coefsFile, name = "coefs")
        }
        else {
            ## process spline information
            coefs <- do.call("cbind", lapply(selist, function(y) y$coefs))

            ## process SummarizedExperiments
            selist <- lapply(selist, function(y) y$se)
            names(selist) <- names(assay(ggd))
        }
        
        knots <- slot(ggs, "knots")[[1]]
        futile.logger::flog.info("Processing done")

        gg <- GenoGAMList(data = selist, id = ggd@id,
                          family = slot(ggs, "family"),
                          design = design(ggd),
                          sizeFactors = sizeFactors(ggd),
                          factorialDesign = colData(ggd),
                          params = modelParams,
                          settings = settings,
                          coefs = coefs,
                          knots = knots,
                          hdf5 = slot(ggd, "hdf5"))
    }

    ## Dataset not split
    ## -----------------
    
    else {

        ## With HDF5 backend
        ## -----------------
        
        if(slot(ggd, "hdf5")) {
            ## create HDF5 file for coefficients
            ## the dimension of the matrix for coefsfile (betas * tile width)
            nfun <- length(.getVars(design(ggd), type = "covar"))
            nbetas <- dim(slot(ggs, "designMatrix"))[2]
            d <- c(nbetas * nfun, length(getIndex(ggd)))

            ## create Coefs file
            seedFile <- .get_seed(assay(ggd))
            chunk <- c(nbetas * nfun, 1)
            coefsFile <- .createH5DF(seedFile, settings, d, chunk, what = "coefs")

            rr <- rowRanges(ggd)
                              
            ## the dimension of the matrix for given chromosome
            d <- c(length(rr), nfun)
            ## create datasets
            h5file <- .createH5DF(seedFile, settings, d, chunk = c(getChunkSize(ggd), nfun))
        }

        ## no HDF5 backend
        ## -----------------
        
        else {
            h5file <- NULL
            coefsFile <- NULL
        }

        qdir <- .init_Queue(h5file)
        res <- BiocParallel::bplapply(ids, .fitGenoGAM, 
                                      data = ggd, init = ggs, coords = coords,
                                      relativeChunks = relativeChunks, h5file = h5file,
                                      chunks = chunks, coefsFile = coefsFile,
                                      qdir = qdir)
        ## remove temporary queue folder and close files
        .end_Queue(qdir)
        
        futile.logger::flog.info("Done")

        ## With HDF5 backend
        ## -----------------
        
        if(slot(ggd, "hdf5")) {
            ## TODO: In Hdf5 write function need to attach log entry to HDF5 dumpLog
            h5fits <- HDF5Array::HDF5Array(h5file, name = "fits")
            h5ses <- HDF5Array::HDF5Array(h5file, name = "ses")

            ## make names, as HDF5 does not store them
            colDataNames <- .makeNames(design(ggd))
            df <- DataFrame(matrix(,length(colDataNames),0))
            rownames(df) <- colDataNames
            
            se <- SummarizedExperiment::SummarizedExperiment(rowRanges = rr,
                                                             assays = list(fits = h5fits,
                                                                           se = h5ses))
            colData(se) <- df
            coefs <- HDF5Array::HDF5Array(coefsFile, name = "coefs")
        }

        ## No HDF5 backend
        ## -----------------

        else {
            futile.logger::flog.info("Processing fits")
            ## build GenoGAM object
            combinedFits <- .transformResults(res, relativeChunks, what = "fits")
            combinedSEs <- .transformResults(res, relativeChunks, what = "se")

            se <- SummarizedExperiment::SummarizedExperiment(rowRanges = SummarizedExperiment::rowRanges(ggd),
                                                             assays = list(fits = combinedFits,
                                                                           se = combinedSEs))
            ## process spline information
            coefs <- sapply(res, function(y) {
                slot(y, "beta")
            })
        }

        knots <- slot(ggs, "knots")[[1]]
        futile.logger::flog.info("Processing done")
        
        gg <- GenoGAM(se, family = slot(ggs, "family"),
                      design = design(ggd),
                      sizeFactors = sizeFactors(ggd),
                      factorialDesign = colData(ggd),
                      params = modelParams,
                      settings = settings,
                      coefs = coefs,
                      knots = knots,
                      hdf5 = slot(ggd, "hdf5"))
    }

    futile.logger::flog.info("Finished")
    return(gg)
}

.fitGenoGAM <- function(id, data, init, coords, relativeChunks = NULL, chunks = NULL, h5file = NULL,
                       coefsFile = NULL, qdir = NULL) {
    suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))

    setup <- .initiate(data, init, coords, id)
    betas <- .estimateParams(setup)

    futile.logger::flog.debug(paste("Beta estimation for tile", id, "took", betas$iterations, "iterations"))
    if(betas$converged == FALSE) {
        futile.logger::flog.warn("Beta estimation did not converge. Increasing the 'maxiter' or 'eps' parameter in 'estimControl' slot in the settings might help, but should be done at own risk.")
    }
    
    slot(setup, "beta") <- betas$par
    slot(setup, "fits") <- .getFits(setup)
    
    slot(setup, "se") <- .compute_SE(setup)
    slot(setup, "params")$id <- id
    
    ## procedure to gradually write results to HDF5, to safe memory footprint
    if(slot(data, "hdf5")) {
        if(any(is.null(chunks), is.null(h5file))) {
            ## An error rather for the developer
            stop("Chunk coordinates missing in the fitting function")
        }

        ## extract fits and SEs
        combinedFits <- .transformResults(list(setup), relativeChunks, what = "fits")
        combinedSEs <- .transformResults(list(setup), relativeChunks, what = "se")

        if(is(data, "GenoGAMDataSetList")) {
            ## normalize ID chromosome-wise, as it is usually based on the entire genome
            chrom <- as.character(GenomeInfoDb::seqnames(getIndex(data)[id,]))
            subindx <- getIndex(data)[GenomeInfoDb::seqnames(getIndex(data)) == chrom,]
            normID <- which(subindx == getIndex(data)[id,])
        }
        else {
            normID <- id
        }

        ## queue write process
        qid <- .queue(qdir)
        ## wait till it's your turn
        
        while(list.files(qdir)[1] != qid){
            Sys.sleep(0.1)
        }

        ## write data
        ## Fits
        rhdf5::h5write(as.matrix(combinedFits), file = h5file, name = "/fits",
                       index = list(start(chunks)[normID]:end(chunks)[normID], 1:ncol(combinedFits)))
        rhdf5::h5write(as.matrix(combinedSEs), file = h5file, name = "/ses",
                       index = list(start(chunks)[normID]:end(chunks)[normID], 1:ncol(combinedFits)))
        ## Coefs
        rhdf5::h5write(betas$par, file = coefsFile, name = "coefs",
                       index = list(1:length(betas$par), id))

        ## Unqueue file by removing
        .unqueue(qid, qdir)
        
        ## reset not needed slots
        slot(setup, "designMatrix") <- new("dgCMatrix")
        slot(setup, "penaltyMatrix") <- new("dgCMatrix")
        slot(setup, "response") <- numeric()
        slot(setup, "offset") <- numeric()
        
        return(NULL)
    }

    return(setup)
}

##########################

## Builds the response vector from GenoGAMDataSet and the given row coordinates
.buildResponseVector <- function(ggd, coords, id) {
    
    if(length(ggd) == 0 | length(coords) == 0) {
        return(numeric())
    }
    
    tile <- coords[id,]
    y <- .subsetByCoords(x = assay(ggd), i = start(tile):end(tile))
    
    Y <- unname(unlist(as.data.frame(y)))
    return(as.integer(Y))
}

## initiates GenoGAMSetup with tile specific data
.initiate <- function(ggd, setup, coords, id) {

    ## build X from X0 and design
    des <- .getDesignFromFormula(design(ggd), colData(ggd))
    X <- as(.blockMatrixFromDesignMatrix(slot(setup, "designMatrix"), des), "dgCMatrix")
    slot(setup, "designMatrix") <- X

    ## build S from S0 and design
    Sdes <- diag(ncol(des))
    S <- as(.blockMatrixFromDesignMatrix(slot(setup, "penaltyMatrix"), Sdes), "dgCMatrix")
    slot(setup, "penaltyMatrix") <- S

    ## initiate response vector and betas
    slot(setup, "response") <- .buildResponseVector(ggd, coords, id)
    numBetas <- dim(slot(setup, "designMatrix"))[2]
    
    if(length(slot(setup, "response")) != 0) {
        ## setup betas as all 1 in normal space = 0s in log space
        betas <- rep(1, numBetas)
        ## set betas to the runmedian of the response
        ## means <- IRanges::runmed(slot(setup, "response"), 11, endrule = "median")
        ## ## take 250 values at equidistant positions
        ## ks <- as.integer(seq(1, length(means), length.out = numBetas))
        ## ## set all zero values to 1, because of log
        ## betas <- means[ks]
        ## betas[which(betas == 0)] <- 1
        slot(setup, "beta") <- matrix(log(betas), numBetas, 1)
    }
        
    ## initialize lambda and theta
    params <- slot(setup, "params")
    if(is.null(params$lambda)) {
        params$lambda <- numBetas
    }
    if(is.null(params$theta)) {
        params$theta <- 1
    }
    slot(setup, "params") <- params

    return(setup)
}

## compute fits from design matrix and estimated betas
.getFits <- function(setup, log = TRUE) {
    if(all(dim(slot(setup, "designMatrix")) == c(0, 0)) |
       all(is.na(slot(setup, "beta")))) {
        return(list())
    }

    ## get design matrix and beta vector
    X <- slot(setup, "designMatrix")
    betas <- slot(setup, "beta")

    ## get row and column indeces of the design matrix
    ## for the respective fits
    des <- slot(setup, "design")
    nSplines <- ncol(des)
    nSamples <- nrow(des)
    dims <- dim(X)
    rows <- list(1:dims[1])
    cols <- list(1:dims[2])

    if(nSplines > 1) {
        cols <- split(1:dims[2], cut(1:dims[2], nSplines, labels = 1:nSplines))
    }
    if(nSamples > 1) {
        rows <- split(1:dims[1], cut(1:dims[1], nSamples, labels = 1:nSamples))
    }

    ## compute the fits by column (splines) and row (samples)
    ## of the design matrix. All fits of the same spline
    ## are combined to a vector to have the same structure as
    ## the response. Each spline is one element of the list, that
    ## gets returned as the result.
    fits <- lapply(1:nSplines, function(j) {
        res <- sapply(1:nSamples, function(i) {
            Xsub <- X[rows[[i]], cols[[j]]]
            beta <- betas[cols[[j]], 1]
            fit <- as.vector(Xsub %*% beta)
            if(!log) {
                fit <- exp(fit)
            }
            return(fit)
        })
        res <- as.vector(res)
        return(res)
    })
    varNames <- .makeNames(slot(setup, "formula"))
    names(fits) <- varNames

    return(fits)
}

## transform the result list of GenoGAMSetup object into a DataFrame
## of either fits or standard errors
.transformResults <- function(x, relativeChunks, id = NULL, what = c("fits", "se")) {
    if(length(relativeChunks) == 0) {
        return(data.table::data.table())
    }
    aslist <- TRUE

    if(is.null(id)) {
        id <- data.table::data.table(index = 1:length(x), listid = 1)
        aslist <- FALSE
        identifiers <- unique(id$listid)
    }

    if(aslist) {
        identifiers <- unique(id$listid)
        res <- vector("list", length(identifiers))
    }
        
    for(ii in identifiers) {
        elements <- id[id$listid == ii,]$index
        
        allData <- lapply(x[elements], function(y) {
            if(length(slot(y, what)) == 0) {
                return(data.table::data.table())
            }
            
            s <- start(relativeChunks[slot(y, "params")$id])
            e <- end(relativeChunks[slot(y, "params")$id])

            ## go by column and subset by colData column
            l <- lapply(1:length(slot(y, what)), function(z) {
                
                des <- slot(y, "design")
                fits <- slot(y, "fits")
                if(length(des) == 0 | length(fits) == 0) {
                    res <- list()
                }
                else {
                    tileLength <- length(fits[[z]])/nrow(des)
                    ## expand colData columns to full length and coerce to boolean
                    desVec <- matrix(as.logical(rep(des, each = tileLength)), ncol = ncol(des))
                    ## take the z list element (= column) and
                    ## subset by the z colData column
                    res <- slot(y, what)[[z]][desVec[,z]]
                    res <- res[s:e]
                }
                return(res)
            })
            names(l) <- names(slot(y, what))
            return(data.table::as.data.table(l))
        })
        combinedData <- data.table::rbindlist(allData)
        varNames <- colnames(combinedData)
        combinedData <- S4Vectors::DataFrame(combinedData)
        colnames(combinedData) <- varNames
        if(aslist) {
            res[[ii]] <- combinedData
        }
        else {
            res <- combinedData
        }
    }

    return(res)
}
