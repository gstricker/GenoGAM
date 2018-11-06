#' general function to compute parameters of a background distribution
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.computeBackground <- function(fit, se, smooth, is_split = FALSE,
                               is_hdf5 = FALSE) {
    
    if(is_split){
        if(is_hdf5) {
            res <- .computeBackground_split_hdf5(fit, se, smooth)
        }
        else {
            res <- .computeBackground_split(fit, se, smooth)
        }
    }
    else {
        if(is_hdf5) {
            res <- .computeBackground_hdf5(fit, se, smooth)
        }
        else {
            res <- .computeBackground_default(fit, se, smooth)
        }
    }

    return(res)
}

#' compute parameters of a background distribution for default GenoGAM object
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.computeBackground_default <- function(fit, se, smooth) {
    ## initializing the result matrix for smooths and mu0/var0 
    res <- data.frame(matrix(,length(smooth), 2), row.names = smooth)
    colnames(res) <- c("mu0", "var0")
        
    for(name in smooth) {
        all <- fit[,name]
        se_all <- se[,name]
        
        ## if genefilter package not present use median
        if(!requireNamespace("genefilter", quietly = TRUE)) {
            futile.logger::flog.info("genefilter package not found, median will be used instead of shorth.")
            mu0 <- median(all, na.rm=TRUE)
        }
        else {
            mu0 <- genefilter::shorth(all, na.rm=TRUE)
        }
        
        mu0_msg <- paste0("mu0 for ", name, " was computed as ", mu0)
        futile.logger::flog.debug(mu0_msg)
        
        left <- all[all <= mu0]
        right <- abs(left - mu0) + mu0
        new_data <- c(left, right)
        var0 <- mad(new_data, na.rm = TRUE)^2
        
        ## assign afterwards for readability
        res[name, 'mu0'] <- mu0
        res[name, 'var0'] <- var0

        var0_msg <- paste0("var0 for ", name, " was computed as ", var0)
        futile.logger::flog.debug(var0_msg)
    }
    return(res)
}

#' compute parameters of a background distribution for non-split HDF5 backend data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.computeBackground_hdf5 <- function(fit, se, smooth) {
    ## initializing the result matrix for smooths and mu0/var0 
    res <- data.frame(matrix(,length(smooth), 2), row.names = smooth)
    colnames(res) <- c("mu0", "var0")

    for(name in smooth) {
        all <- fit[,name, drop = FALSE]
        se_all <- se[,name, drop = FALSE]
        
        ## if genefilter package not present use median
        if(!requireNamespace("genefilter", quietly = TRUE)) {
            futile.logger::flog.info("genefilter package not found, median will be used instead of shorth.")
            mu0 <- .block_APPLY(all, .median_hdf5())
        }
        else {
            mu0 <- .block_APPLY(all, .shorth_hdf5())
        }
        
        mu0_msg <- paste0("mu0 for ", name, " was computed as ", mu0)
        futile.logger::flog.debug(mu0_msg)

        var0 <- .block_APPLY(all, .var0_hdf5(), mu0)

        ## assign afterwards for readability
        res[name, 'mu0'] <- mu0
        res[name, 'var0'] <- var0

        var0_msg <- paste0("var0 for ", name, " was computed as ", var0)
        futile.logger::flog.debug(var0_msg)
    }
    return(res)
}

#' compute parameters of a background distribution for split data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.computeBackground_split <- function(fit, se, smooth) {
    ## initializing the result matrix for smooths and mu0/var0 
    res <- data.frame(matrix(,length(smooth), 2), row.names = smooth)
    colnames(res) <- c("mu0", "var0")

    for(name in smooth) {
        all <- lapply(fit, function(x) x[,name])
        se_all <- lapply(se, function(x) x[,name])
        
        ## if genefilter package not present use median
        if(!requireNamespace("genefilter", quietly = TRUE)) {
            futile.logger::flog.info("genefilter package not found, median will be used instead of shorth.")
            mu0 <- .mu0_split(all, type = "median")
        }
        else {
            mu0 <- .mu0_split(all, type = "shorth")
        }
        
        mu0_msg <- paste0("mu0 for ", name, " was computed as ", mu0)
        futile.logger::flog.debug(mu0_msg)

        var0 <- .var0_split(all, mu0)

        ## assign afterwards for readability
        res[name, 'mu0'] <- mu0
        res[name, 'var0'] <- var0

        var0_msg <- paste0("var0 for ", name, " was computed as ", var0)
        futile.logger::flog.debug(var0_msg)
    }
    return(res)
}

#' compute parameters of a background distribution for split HDF5 backend data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.computeBackground_split_hdf5 <- function(fit, se, smooth) {
    ## initializing the result matrix for smooths and mu0/var0 
    res <- data.frame(matrix(,length(smooth), 2), row.names = smooth)
    colnames(res) <- c("mu0", "var0")

    for(name in smooth) {
        all <- lapply(fit, function(x) x[,name, drop = FALSE])
        se_all <- lapply(se, function(x) x[,name, drop = FALSE])
        
        ## if genefilter package not present use median
        if(!requireNamespace("genefilter", quietly = TRUE)) {
            futile.logger::flog.info("genefilter package not found, median will be used instead of shorth.")
            mu0 <- .mu0_split_hdf5(all, type = "median")
        }
        else {
            mu0 <- .mu0_split_hdf5(all, type = "shorth")
        }
        
        mu0_msg <- paste0("mu0 for ", name, " was computed as ", mu0)
        futile.logger::flog.debug(mu0_msg)

        var0 <- .var0_split_hdf5 (all, mu0)

        ## assign afterwards for readability
        res[name, 'mu0'] <- mu0
        res[name, 'var0'] <- var0

        var0_msg <- paste0("var0 for ", name, " was computed as ", var0)
        futile.logger::flog.debug(var0_msg)
    }
    return(res)
}

## HDF5 helper functions (mostly for block_APPLY)
## get median of medians
## chromosome dependent medians vary too much, especially for small organisms
## ===============================================

#' Compute mu0 for split data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.mu0_split <- function(x, type = c("shorth", "median")) {
    ## with median/shorth, adjustment for chromosome length
    ## has to be done by replicating the medians of chromosomes
    ## proportionally to their length. The following is a way to determine
    ## how often to replicate each chromosome median/shorth with the smallest
    ## chromosome being replicated only once
    type = match.arg(type)
    lengths <- sapply(x, length)
    rep_per_chr <- round(lengths/min(lengths))

    ## initializing final object
    idx_start <- c(1, cumsum(rep_per_chr)[-length(x)] + 1)
    idx_end <- cumsum(rep_per_chr)
    sqn <- numeric(sum(rep_per_chr))
    
    for(ii in 1:length(x)) {
        idx <- idx_start[ii]:idx_end[ii]
        if(type == "shorth") {
            sqn[idx] <- rep(genefilter::shorth(x[[ii]], na.rm = TRUE),
                            rep_per_chr[ii])
        }
        else {
            sqn[idx] <- rep(median(x[[ii]], na.rm = TRUE),
                            rep_per_chr[ii])
        }
    }

    ## median/shorth of medians/shorths
    if(type == "shorth") {
        res <- tryCatch({
            genefilter::shorth(sqn, na.rm=TRUE)       
        }, error = function(e) {
            median(sqn, na.rm = TRUE)
        })
    }
    else {
        res <- median(sqn, na.rm = TRUE)
    }
    return(res)
}

#' Compute mu0 for split and HDf5 backend data
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.mu0_split_hdf5 <- function(x, type = c("shorth", "median")) {
    ## with median/shorth, adjustment for chromosome length
    ## has to be done by replicating the medians of chromosomes
    ## proportionally to their length. The following is a way to determine
    ## how often to replicate each chromosome median/shorth with the smallest
    ## chromosome being replicated only once
    type = match.arg(type)
    lengths <- sapply(x, length)
    rep_per_chr <- round(lengths/min(lengths))

    ## initializing final object
    idx_start <- c(1, cumsum(rep_per_chr)[-length(x)] + 1)
    idx_end <- cumsum(rep_per_chr)
    sqn <- numeric(sum(rep_per_chr))
    
    for(ii in 1:length(x)) {
        idx <- idx_start[ii]:idx_end[ii]
        if(type == "shorth") {
            sqn[idx] <- rep(.block_APPLY(x[[ii]], .shorth_hdf5()),
                            rep_per_chr[ii])
        }
        else {
            sqn[idx] <- rep(.block_APPLY(x[[ii]], .median_hdf5()),
                            rep_per_chr[ii])
        }
    }

    ## median/shorth of medians/shorths
    if(type == "shorth") {
        res <- tryCatch({
            genefilter::shorth(sqn, na.rm=TRUE)       
        }, error = function(e) {
            median(sqn, na.rm = TRUE)
        })
    }
    else {
        res <- median(sqn, na.rm = TRUE)
    }
    return(res)
}

#' Compute var0 for split data by weighted mean of MADs
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.var0_split <- function(x, mu0) {
    lengths <- sapply(x, length)
    var0_vec <- numeric(length(x))
    for(ii in 1:length(x)) {
        left <- x[[ii]][x[[ii]] <= mu0]
        right <- abs(left - mu0) + mu0
        new_data <- c(left, right)
        ## the last part of the equation is the weight of a var0 to adjust
        ## for chromosome length
        var0_vec[ii] <- mad(new_data, na.rm = TRUE)^2 * (lengths[ii]/sum(lengths))
    }
    ## taking the weighted mean (weights already in)
    var0 <- sum(var0_vec)
    return(var0)
}

#' Compute var0 for split and HDF5 backend data through block_APPLY
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.var0_split_hdf5 <- function(x, mu0) {
    lengths <- sapply(x, length)
    var0_vec <- numeric(length(x))
    for(ii in 1:length(x)) {
        ## the last part of the equation is the weight of a var0 to adjust
        ## for chromosome length
        var0_chr <- .block_APPLY(x[[ii]], .var0_hdf5(), mu0)
        var0_vec[ii] <- var0_chr * (lengths[ii]/sum(lengths))
    }
    ## taking the weighted mean (weights already in)
    var0 <- sum(var0_vec)
    return(var0)
}

#' The APPLY-COMBINE function for var0 computation on a single range.
#' Used in higher functions to compute complete var0
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.var0_hdf5 <- function(block, mu0) {
    
    APPLY <- function(block, mu0) {
        tmp <- as.vector(block, mode = "numeric")
        left <- block[block <= mu0]
        right <- abs(left - mu0) + mu0
        new_data <- c(left, right)
        ## the last part of the equation is the weight of a var0 to adjust
        ## for chromosome length
        res <- mad(new_data, na.rm = TRUE)^2
        return(res)
    }

    COMBINE <- function(x) {
        median(x, na.rm = TRUE)
    }

    return(list(APPLY = APPLY, COMBINE = COMBINE))
}

#' The APPLY-COMBINE function for median computation on a single range.
#' Used in higher functions to compute complete median
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.median_hdf5 <- function(block) {
    
    APPLY <- function(block) {
        tmp <- as.vector(block, mode = "numeric")
        res <- median(tmp, na.rm = TRUE)
        return(res)
    }

    COMBINE <- function(x) {
        median(x, na.rm = TRUE)
    }

    return(list(APPLY = APPLY, COMBINE = COMBINE))
}

#' The APPLY-COMBINE function for shorth computation on a single range.
#' Used in higher functions to compute complete shorth
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.shorth_hdf5 <- function(block) {
    
    APPLY <- function(block) {
        tmp <- as.vector(block, mode = "numeric")
        res <- tryCatch({
            genefilter::shorth(tmp, na.rm = TRUE)
        }, error = function(e) {
            median(tmp, na.rm = TRUE)
        })
        return(res)
    }

    COMBINE <- function(x) {
        tryCatch({
            genefilter::shorth(x, na.rm = TRUE)
        }, error = function(e) {
            median(tmp, na.rm = TRUE)
        })
    }

    return(list(APPLY = APPLY, COMBINE = COMBINE))
}
