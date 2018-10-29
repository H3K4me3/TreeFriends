## This script will produce a table showing SNP and their mapped
## nucleotides on other genomes.


library(here)

library(Rsamtools)
library(rtracklayer)
library(VariantAnnotation)

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Ptroglodytes.UCSC.panTro5)
library(BSgenome.GGorilla.UCSC.gorGor5)
library(BSgenome.PAbelii.UCSC.ponAbe2)
library(BSgenome.Mmulatta.UCSC.rheMac8)


library(BiocParallel)

setwd(here())

slidingRanges <- local({
    slidingRanges <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
    slidingRanges <- slidingRanges[seqnames(slidingRanges) %in%
        names(genome(VcfFile(here("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz"))))]
    slidingRanges <- slidingWindows(slidingRanges, width = 1000000, step = 1000000)
    unname(unlist(slidingRanges, use.names = FALSE))
})
slidingRanges
#slidingRanges <- slidingRanges[24]
cat(sprintf("======== In total %s sliding ranges ========\n", length(slidingRanges)))
table(seqnames(slidingRanges))[table(seqnames(slidingRanges)) != 0]


read_vcf <- function(range) {
    stopifnot(length(range) == 1)
    
    # Drop unused seqlevels
    seqlevels(range) <- as.character(unique(seqnames(range)))
    
    vcf <- readVcf(
        TabixFile(here("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz")),
        param = ScanVcfParam(which = range)
    )
    
    ans <- rowRanges(vcf)
    
    # Assign strand to positive strand
    stopifnot(all(strand(ans) == "*"))
    strand(ans) <- "+"
    
    # Remove some columns
    ans$QUAL <- NULL
    ans$paramRangeID <- NULL
    ans$FILTER <- NULL
    
    # Remove names
    mcols(ans) <- cbind(DataFrame(id = names(ans)), mcols(ans))
    names(ans) <- NULL
    
    # Only select ones that has start(ans) lays inside
    # the range to avoid duplicated results.
    ans[start(ans) >= start(range) & start(ans) <= end(range)]
    
    # Unlist ALT
    stopifnot(is(ans$ALT, "DNAStringSetList"))
    stopifnot(all(lengths(ans$ALT) == 1))
    ans$ALT <- unlist(ans$ALT)
    
    ans
}


lift_over <- local({
    suppressWarnings({
        chain.panTro5 <- import.chain("raw_data/chainfiles/hg38ToPanTro5.over.chain")
        chain.gorGor5 <- import.chain("raw_data/chainfiles/hg38ToGorGor5.over.chain")
        chain.ponAbe2 <- import.chain("raw_data/chainfiles/hg38ToPonAbe2.over.chain")
        chain.rheMac8 <- import.chain("raw_data/chainfiles/hg38ToRheMac8.over.chain")
    })
    
    function(gr) {
        lift <- function(gr, chain) {
            if (is.character(chain))
                chain <- import.chain(chain)
            
            grl.target <- local({
                mcols(gr) <- NULL
                liftOver(gr, chain)
            })
            grl.target
        }
        queryString <- function(grl, bsgenome) {
            gr <- unlist(unname(grl))
            dna <- DNAStringSet(Views(bsgenome, gr))
            relist(dna, unname(grl))
        }
        
        gr$loc.panTro5 <- lift(gr, chain.panTro5)
        gr$dna.panTro5 <- queryString(gr$loc.panTro5, BSgenome.Ptroglodytes.UCSC.panTro5)
        
        gr$loc.gorGor5 <- lift(gr, chain.gorGor5)
        gr$dna.gorGor5 <- queryString(gr$loc.gorGor5, BSgenome.GGorilla.UCSC.gorGor5)
        
        gr$loc.ponAbe2 <- lift(gr, chain.ponAbe2)
        gr$dna.ponAbe2 <- queryString(gr$loc.ponAbe2, BSgenome.PAbelii.UCSC.ponAbe2)
        
        gr$loc.rheMac8 <- lift(gr, chain.rheMac8)
        gr$dna.rheMac8 <- queryString(gr$loc.rheMac8, BSgenome.Mmulatta.UCSC.rheMac8)
        
        gr
    }
})


library(SNPlocs.Hsapiens.dbSNP151.GRCh38)

make_snp_table <- function(gr) {
    gr <- gr[width(gr) == 1]
    ans <- tibble::tibble(
        Loc = as.character(gr),
        REF = as.character(gr$REF),
        ALT = as.character(gr$ALT)
    )
    stopifnot(all(unique(lengths(gr$dna.panTro5)) %in% c(0, 1)))
    stopifnot(all(unique(lengths(gr$dna.gorGor5)) %in% c(0, 1)))
    stopifnot(all(unique(lengths(gr$dna.ponAbe2)) %in% c(0, 1)))
    stopifnot(all(unique(lengths(gr$dna.rheMac8)) %in% c(0, 1)))
    
    to_char <- function(dnasetl) {
        stopifnot(all(unique(lengths(dnasetl)) %in% c(0, 1)))
        ans <- replicate(length(dnasetl), NA_character_)
        ans[lengths(dnasetl) == 1] <- as.character(unlist(dnasetl))
        ans
    }
    ans$panTro5 <- to_char(gr$dna.panTro5)
    ans$gorGor5 <- to_char(gr$dna.gorGor5)
    ans$ponAbe2 <- to_char(gr$dna.ponAbe2)
    ans$rheMac8 <- to_char(gr$dna.rheMac8)
    
    ans
}

annotate_rsid <- function(table) {
    loc <- GPos(table$Loc)
    seqlevelsStyle(loc) <- seqlevelsStyle(SNPlocs.Hsapiens.dbSNP151.GRCh38)[1]
    stopifnot(length(unique(seqnames(loc))) == 1)
    snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP151.GRCh38, range(loc))
    
    rsid <- vector("list", length(loc))
    hits <- findMatches(ranges(loc), ranges(snps))
    
    q_hits <- queryHits(hits)
    tmp <- split(snps$RefSNP_id[subjectHits(hits)], q_hits)
    rsid[as.integer(names(tmp))] <- tmp
    table$rsid <- rsid
    table
}


write_table <- function(x, path) {
    x$rsid <- sapply(x$rsid, function(x) {
        if (length(x) == 1)
            x[1]
        else
            paste(x, collapse = ",")
    })
    readr::write_tsv(x, path, na = ".")
}


res_dir <- here("results", paste0("snptable_", format(Sys.time(), "%Y-%m-%d_%H:%M")))
if (!dir.exists(res_dir))
    dir.create(res_dir, recursive = TRUE)

bplapply(seq_along(slidingRanges), BPPARAM = MulticoreParam(),
    function(i) {
        d <- slidingRanges[i]
        d <- read_vcf(d)
        d <- lift_over(d)
        d <- make_snp_table(d)
        d <- annotate_rsid(d)
        path <- file.path(res_dir, sprintf("%04d_%s.txt.gz", i,
                                           format(Sys.time(), "%Y-%m-%d_%H:%M")))
        write_table(d, path)
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



