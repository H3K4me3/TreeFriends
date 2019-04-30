
library(dplyr)
library(here)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(GenomicRanges)

stopifnot(
    tools::md5sum(here("raw_data/Site_frequency_table_apr25.txt")) == "0a8424afd597252c61ee0e2522dbd4e6"
)

frequency_table <- readr::read_tsv(here("raw_data/Site_frequency_table_apr25.txt"))

## 23 is chrX
## 24 is chrY
## 25 is chrX
## 26 is chrM
frequency_table <- dplyr::mutate(frequency_table, CHR = paste0("chr", CHR)) %>%
    dplyr::mutate(CHR = if_else(CHR == "chr23", "chrX", CHR)) %>%
    dplyr::mutate(CHR = if_else(CHR == "chr24", "chrY", CHR)) %>%
    dplyr::mutate(CHR = if_else(CHR == "chr25", "chrX", CHR)) %>%
    dplyr::mutate(CHR = if_else(CHR == "chr26", "chrM", CHR))
frequency_table

loc <- GRanges(frequency_table$CHR, IRanges(frequency_table$pos_hg38, width = 1))
loc <- as.data.frame(loc)
loc <- loc %>% dplyr::mutate(seqnames = as.character(seqnames)) %>%
    dplyr::select(-strand)
loc


search_db <- function(seqnames, location) {
    stopifnot(is.character(seqnames))
    location <- as.integer(location)
    db <- src_sqlite(path = here("results/snpdb.sqlite3"))
    snp <- tbl(db, "snp")
    sel <- dplyr::select(
        snp, seqnames, start, `hg-panTro`, `hg-panTro-gorGor` = `hg-pan-gor`,
        `hg-panTro-gorGor-ponAbe` = `ponAbe-else`,
        `hg-panTro-gorGor-ponAbe-rheMac` = root,
        hg = hg, panTro = panTro, gorGor = gorGor, ponAbe = ponAbe, rheMac = rheMac
    )
    
    key <- tibble::tibble(seqnames = seqnames, start = location)
    res <- left_join(key, sel, copy = TRUE)
    res
}

ans <- search_db(loc$seqnames, loc$start)

readr::write_tsv(ans, here("results/ExtractAncestralAllele_apr25.txt"))

# source(here("lib/liftover.R"))
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome.Ptroglodytes.UCSC.panTro5)
# library(BSgenome.GGorilla.UCSC.gorGor5)
# library(BSgenome.PAbelii.UCSC.ponAbe2)
# library(BSgenome.Mmulatta.UCSC.rheMac8)
# 
# suppressWarnings({
#     chain.panTro5 <- import.chain(here("raw_data/chainfiles/hg38ToPanTro5.over.chain"))
#     chain.gorGor5 <- import.chain(here("raw_data/chainfiles/hg38ToGorGor5.over.chain"))
#     chain.ponAbe2 <- import.chain(here("raw_data/chainfiles/hg38ToPonAbe2.over.chain"))
#     chain.rheMac8 <- import.chain(here("raw_data/chainfiles/hg38ToRheMac8.over.chain"))
# })
# 
# d1 <- (lift_over(loc, BSgenome.Ptroglodytes.UCSC.panTro5, chain.panTro5))
# d2 <- (lift_over(loc, BSgenome.GGorilla.UCSC.gorGor5, chain.gorGor5))
# d3 <- (lift_over(loc, BSgenome.PAbelii.UCSC.ponAbe2, chain.ponAbe2))
# d4 <- (lift_over(loc, BSgenome.Mmulatta.UCSC.rheMac8, chain.rheMac8))
# 
# lift_res <- cbind(d1, d2, d3, d4)[,c(2,4,6,8)]
# 
# rm(d1); rm(d2); rm(d3); rm(d4)
# 
# lift_res <- cbind(loc, lift_res)
# lift_res

