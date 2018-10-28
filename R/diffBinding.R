#' Compute significance for given regions
#'
#' For a given set of regions, region-wise pvalues and FDR is computed
#'
#' @param fit A GenoGAM object containing the fit
#' @param regions A GRanges object of regions of interest
#' @param smooth Which fit should be used. The names should be equivalent to the column names
#' of the object. Lookup with \code{colnames(my_GenoGAM_object)}
#' @return The GRanges object from the 'region' parameter extended by two
#' columns: pvalue and FDR
#' @details For a given set of regions, region-wise pvalues are computed by applying
#' familywise hochberg correction and taking the minimal p-value. FDR is computed by further
#' applying Benjamini-Hochberg correction.
#' @examples
#' ## make test GenoGAM
#' gg <- makeTestGenoGAM()
#' ## make region
#' region <- GRanges("chrXYZ", IRanges(c(2000, 4000, 6000), c(3000, 5000, 9000)))
#' res <- computeRegionSignificance(gg, region)
#' res
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
computeRegionSignificance <- function(fit, regions, smooth = NULL) { 
    futile.logger::flog.info("Estimating region p-values and FDR")

    ## which GenoGAM backend is present
    if(is(fit, "GenoGAMList")) {
        is_split <- TRUE
    }
    else {
        if(is(fit, "GenoGAM")) {
            is_split <- FALSE
        }

        else {
            stop("Wrong class submitted")
        }
    }

    ## if smooth unknown do for all smooths
    if(is.null(smooth)) {
        smooth <- colnames(fit)
    }

    ## first subset by regions
    fit_regions <- fit[regions]
    ## compute pvalues on reduced fit
    computeSignificance(fit_regions, log.p = TRUE)

    ## no need for HDF5 conditions as it does not affect the computation
    if(is_split){
        pvals <- .get_pvals_split(fit_regions, regions)
    }
    else {
        pvals <- .get_pvals_default(fit_regions, regions)
    }

    futile.logger::flog.info("Performing multiple correction")
    res <- vector("list", length(smooth))
    names(res) <- smooth
    for(name in smooth) {
        res[[name]] <- regions
        ## use setnames to get around the problem that subsetting data.table
        ## by quoted column names which are within a function is complicated
        data.table::setnames(pvals, name, "p")
        pv = pvals[, min(p.adjust(exp(p), method="hochberg")), by = region]
        data.table::setnames(pvals, "p", name)

        res[[name]]$pvalue = NA
        res[[name]]$pvalue[pv$region] = pv$V1
        res[[name]]$FDR = p.adjust(res[[name]]$pvalue, method="BH")
    }

    futile.logger::flog.info("DONE")
    return(res)
}

.get_pvals_default <- function(fit, regions) {
    ## Determine indices in DataFrame. This is only necessary since
    ## some regions might overlap
    regions_group <- IRanges::findOverlaps(rowRanges(fit), regions)
    gg <- fit[S4Vectors::queryHits(regions_group)]
    
    res <- data.table::data.table(data.frame(pvalue(gg)))
    ## during conversion the names get changed, so change back
    data.table::setnames(res, names(res), colnames(gg))
    res$region <- as.factor(S4Vectors::subjectHits(regions_group))
    
    return(res)
}

.get_pvals_split <- function(fit, regions) {
        
    ## Determine indices in each DataFrame. This is only necessary since
    ## some regions might overlap
    regions_group <- lapply(rowRanges(fit), function(x) {
        IRanges::findOverlaps(x, regions)
    })

    ## Subset pvalue DataFrame according to indices 
    gg <- lapply(names(regions_group), function(chr) {
        pvalue(fit)[[chr]][S4Vectors::queryHits(regions_group[[chr]]),]        
    })

    gg <- do.call("rbind", gg)
    res <- data.table::data.table(data.frame(gg))
    ## during conversion the names get changed, so change back
    data.table::setnames(res, names(res), colnames(gg))
    ## flatten regions_group to add as column
    regions_group <- unlist(unname(sapply(regions_group, S4Vectors::subjectHits)))
    res$region <- as.factor(regions_group)

    return(res)
}
