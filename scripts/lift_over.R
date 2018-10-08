
library(here)

library(Rsamtools)
library(rtracklayer)
library(VariantAnnotation)

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Ptroglodytes.UCSC.panTro5)
library(BSgenome.GGorilla.UCSC.gorGor5)
library(BSgenome.PAbelii.UCSC.ponAbe2)

# In order to get snp id
#library(SNPlocs.Hsapiens.dbSNP151.GRCh38)

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
cat(sprintf("======== In total %s sliding ranges ========\\n", length(slidingRanges)))
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
    ans$QUAL <- NULL
    ans$paramRangeID <- NULL
    stopifnot(all(strand(ans) == "*"))
    strand(ans) <- "+"
    ans
    
    # Only select ones that has start(ans) lays inside
    # the range to avoid duplicated results.
    ans[start(ans) >= start(range) & start(ans) <= end(range)]
}

lift_over <- local({
    suppressWarnings({
        chain.panTro5 <- import.chain("raw_data/chainfiles/hg38ToPanTro5.over.chain")
        chain.gorGor5 <- import.chain("raw_data/chainfiles/hg38ToGorGor5.over.chain")
        chain.ponAbe2 <- import.chain("raw_data/chainfiles/hg38ToPonAbe2.over.chain")
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
        
        gr
    }
})

# TODO: get snp id

if (!dir.exists(here("lift_over_results")))
    dir.create(here("lift_over_results"))

bplapply(seq_along(slidingRanges), BPPARAM = MulticoreParam(),
    function(i) {
        gr <- slidingRanges[i]
        gr <- read_vcf(gr)
        gr <- lift_over(gr)
        saveRDS(gr, here(sprintf("lift_over_results/%04d_%s.rds", i,
                             format(Sys.time(), "%Y-%m-%d_%H:%M"))))
        return(TRUE)
    }
)



