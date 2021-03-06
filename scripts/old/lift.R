
library(here)

library(Rsamtools)
library(rtracklayer)
library(VariantAnnotation)

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Ptroglodytes.UCSC.panTro5)
library(BSgenome.GGorilla.UCSC.gorGor5)
library(BSgenome.PAbelii.UCSC.ponAbe2)
library(BSgenome.Mmulatta.UCSC.rheMac8)

#library(BiocParallel)

slidingGenome <- function(width = 1000000, chr = NULL) {
    vcffile <- VcfFile(here("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz"))
    # Get chromosome and chromosome length
    slidingRanges <- GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))
    # Only keep necessary chromosomes
    slidingRanges <- slidingRanges[seqnames(slidingRanges) %in% names(genome(vcffile))]
    
    slidingRanges <- slidingWindows(slidingRanges, width = width, step = width)
    
    # Restrict to given chromosomes
    if (!is.null(chr))
        slidingRanges <- slidingRanges[seqnames(slidingRanges) %in% chr]
    
    unname(unlist(slidingRanges, use.names = FALSE))
}


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
    ans <- ans[start(ans) >= start(range) & start(ans) <= end(range)]
    
    # Unlist ALT
    stopifnot(is(ans$ALT, "DNAStringSetList"))
    stopifnot(all(lengths(ans$ALT) == 1))
    ans$ALT <- unlist(ans$ALT)
    
    ans
}


lift_over <- local({
    suppressWarnings({
        chain.panTro5 <- import.chain(here("raw_data/chainfiles/hg38ToPanTro5.over.chain"))
        chain.gorGor5 <- import.chain(here("raw_data/chainfiles/hg38ToGorGor5.over.chain"))
        chain.ponAbe2 <- import.chain(here("raw_data/chainfiles/hg38ToPonAbe2.over.chain"))
        chain.rheMac8 <- import.chain(here("raw_data/chainfiles/hg38ToRheMac8.over.chain"))
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


make_snp_table <- function(gr) {
    ## Filter the insertion and deletion, only keep SNPs
    gr <- gr[width(gr$REF) == 1 & width(gr$ALT) == 1]
    
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

#library(SNPlocs.Hsapiens.dbSNP151.GRCh38)

#annotate_rsid <- function(table) {
#    loc <- GPos(table$Loc)
#    seqlevelsStyle(loc) <- seqlevelsStyle(SNPlocs.Hsapiens.dbSNP151.GRCh38)[1]
#    stopifnot(length(unique(seqnames(loc))) == 1)
#    snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP151.GRCh38, range(loc))
#    
#    rsid <- vector("list", length(loc))
#    hits <- findMatches(ranges(loc), ranges(snps))
#    
#    q_hits <- queryHits(hits)
#    tmp <- split(snps$RefSNP_id[subjectHits(hits)], q_hits)
#    rsid[as.integer(names(tmp))] <- tmp
#    table$rsid <- rsid
#    table
#}



