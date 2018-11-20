#!/usr/bin/env Rscript

## This script will produce a table showing SNP and their mapped
## nucleotides on other genomes.

library(here)
setwd(here())

source(here('lib/lift.R'))

slidingRanges <- local({
    slidingRanges <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
    slidingRanges <- slidingRanges[seqnames(slidingRanges) %in%
        names(genome(VcfFile(here("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz"))))]
    slidingRanges <- slidingWindows(slidingRanges, width = 1000000, step = 1000000)
    unname(unlist(slidingRanges, use.names = FALSE))
})
#slidingRanges <- slidingRanges[24]
cat(sprintf("======== In total %s sliding ranges ========\n", length(slidingRanges)))
table(seqnames(slidingRanges))[table(seqnames(slidingRanges)) != 0]



write_table <- function(x, path) {
    #x$rsid <- sapply(x$rsid, function(x) {
    #    if (length(x) == 1)
    #        x[1]
    #    else
    #        paste(x, collapse = ",")
    #})
    readr::write_tsv(x, path, na = ".")
}


res_dir <- here("results", paste0("snptable_", format(Sys.time(), "%Y-%m-%d_%H:%M")))
if (!dir.exists(res_dir))
    dir.create(res_dir, recursive = TRUE)

cat(("======== Start processing the data  ========\n"))

library(BiocParallel)

bplapply(seq_along(slidingRanges), BPPARAM = MulticoreParam(),
    function(i) {
        cat(sprintf("-- process %s --", i))
        d <- slidingRanges[i]
        d <- read_vcf(d)
        d <- lift_over(d)
        d <- make_snp_table(d)
        #d <- annotate_rsid(d)
        path <- file.path(res_dir, sprintf("%04d_%s.txt.gz", i,
                                           format(Sys.time(), "%Y-%m-%d_%H:%M")))
        write_table(d, path)
        cat(sprintf("-- finish  %s --", i))
    }
)


#if (!dir.exists(here("results/lift_over")))
#    dir.create(here("results/lift_over"), recursive = TRUE)
#bplapply(seq_along(slidingRanges), BPPARAM = MulticoreParam(),
#    function(i) {
#        gr <- slidingRanges[i]
#        gr <- read_vcf(gr)
#        gr <- lift_over(gr)
#        saveRDS(gr, here(sprintf("results/lift_over/%04d_%s.rds", i,
#                             format(Sys.time(), "%Y-%m-%d_%H:%M"))))
#        return(TRUE)
#    }
#)



