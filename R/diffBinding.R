#' Compute significance for given regions
#'
#' For a given set of regions, region-wise pvalues and FDR is computed
#'
#' @param fit A GenoGAM object containing the fit
#' @param regions A GRanges object of regions of interest
#' @return The GRanges object from the 'region' parameter extended by two columns: pvalue and FDR
#' @details For a given set of regions, region-wise pvalues are computed by applying
#' familywise hochberg correction and taking the minimal p-value. FDR is computed by further
#' applying Benjamini-Hochberg correction.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
computeRegionSignificance <- function(fit, regions) {
    futile.logger::flog.info("Estimating region p-values and FDR")
    
    ov <- data.table::data.table(view(fit, ranges = regions))
    regions_group <- IRanges::findOverlaps(rowRanges(fit), regions)
    ov$gene <- as.factor(subjectHits(regions_group))
    
    data.table::setnames(ov, "pvalue.s(x):genotype", "pvalue")
    pv = ov[, min(p.adjust(pvalue, method="hochberg")), by = gene]
    
    regions$pvalue = NA
    regions$pvalue[pv$gene] = pv$V1
    regions$FDR = p.adjust(regions$pvalue, method="BH")
    futile.logger::flog.info("Done")
    return(regions)
}
