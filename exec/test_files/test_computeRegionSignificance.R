#################################
## A test script new functions ##
#################################

## test function: computeRegionSignificance
library(rtracklayer)
gffPath <- "/s/genomes/sacCer/S288C_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.pure.gff"
gff <- import(gffPath)
gff <- gff[gff$type == "gene",]
gff$name <- ifelse(!is.na(gff$gene), gff$gene, gff$Name)
seqlevels(gff)[17] <- c("chrM")
genes <- gff[, c("ID", "name")]

what <- 'genotype'
##function
diffregions <- computeRegionSignificance(fit, genes, what = what)

## manual
gene_group <- findOverlaps(fit@positions, genes)
ov <- data.table::data.table(view(fit, ranges = genes))
ov$gene <- as.factor(subjectHits(gene_group))
## apply familywise hochberg correction to get region-wise pvalues
## for convenience change column name
data.table::setnames(ov, paste("pvalue.s(x)", what, sep = ":"), "pvalue")
genes_pv = ov[, min(p.adjust(pvalue, method="hochberg")), by = gene]
genes$pvalue = NA
genes$pvalue[genes_pv$gene] = genes_pv$V1
## apply BH correction to get FDR per region
genes$corPvalue = p.adjust(genes$pvalue, method="BH")

comp1 <- all.equal(genes$pvalue, diffregions$pvalue)
comp2 <- all.equal(genes$corPvalue, diffregions$FDR)

(if(comp1 & comp2) {
    cat("SUCCESFUL\n")
}
else {
    cat("ERROR\n")
})
    

