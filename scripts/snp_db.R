
library(here)

library(Rsamtools)
library(rtracklayer)
library(VariantAnnotation)

library(Biostrings)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Ptroglodytes.UCSC.panTro5)
library(BSgenome.GGorilla.UCSC.gorGor5)
library(BSgenome.PAbelii.UCSC.ponAbe2)
library(BSgenome.Mmulatta.UCSC.rheMac8)
#library(RSQLite)

VCF_FILE_LOC = here("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz")
DB_LOC = here("results/snpdb.sqlite")

slidingGenome <- function(width = 1000000, chr = NULL) {
    vcffile <- VcfFile(VCF_FILE_LOC)
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

extractSNPFromVcf <- function(vcf, range = NULL) {
    if (!is.null(range)) {
        stopifnot(length(range) == 1)
        # Drop unused seqlevels
        seqlevels(range) <- as.character(unique(seqnames(range)))
    }
    
    if (is.character(vcf))
        vcf <- VariantAnnotation::VcfFile(vcf)
    if (is(vcf, "VcfFile")) {
        if (is.null(range))
            vcf <- VariantAnnotation::readVcf( vcf, row.names = FALSE )
        else
            vcf <- VariantAnnotation::readVcf( vcf, row.names = FALSE,
                param = ScanVcfParam(which = range)
            )
    }
    stopifnot(is(vcf, "VCF"))
    
    #vcf <- expand(vcf)
    gr <- rowRanges(vcf)
    
    ## Remove row names
    names(gr) <- NULL
    
    ## Add extra columns
    mcols(gr) <- cbind(mcols(gr), info(vcf))
    
    ## Remove some columns that I do not care about
    ## - AF can be calculated from AN and AC
    for (colname in c("paramRangeID", "QUAL", "FILTER", "AF", "VRT"))
        mcols(gr)[[colname]] <- NULL
    
    ## Unlist some columns
    for (colname in c("ALT", "AC", "Het", "Hom")) {
        stopifnot(is(mcols(gr)[[colname]], "List"))
        stopifnot(unique(lengths(mcols(gr)[[colname]])) == 1)
        mcols(gr)[[colname]] <- unlist(mcols(gr)[[colname]])
    }
    
    ## Fix strand
    stopifnot(unique(strand(gr)) == "*")
    strand(gr) <- "+"
    
    # Filter only SNPs
    # FIXME: Note that when we calculate the frequency of the reference allele,
    # we need to taking in consideration of the deletino... It's complex...
    gr <- gr[width(gr) == 1 & width(gr$ALT) == 1]
    
    ## For sanity, ensure all results is within the range
    if (!is.null(range))
        stopifnot(countOverlaps(gr, range, type = "within") == 1)
    ## Otherwise, we::
    # Only select ones that has start(ans) lays inside
    # the range to avoid duplicated results.
    #ans <- ans[start(ans) >= start(range) & start(ans) <= end(range)]
    
    gr
}

#writeSNPToDB <- function(snp, db) {
#    if (is.character(db)) {
#        db <- dbConnect(RSQLite::SQLite(), db)
#        on.exit(dbDisconnect(db))
#    }
#    stopifnot(is(db, "SQLiteConnection"))
#    stopifnot(is(snp, "GRanges"))
#    df <- tibble::as_tibble(as.data.frame(snp))
#    
#    # Drop some columns
#    stopifnot(all(df$width == 1L))
#    df[["width"]] <- NULL
#    
#    if (! "SNP" %in% dbListTables(db)) {
#        message("Creating table 'SNP' in db")
#        dbCreateTable(db, name = "SNP", fields = df)
#    }
#    
#    dbAppendTable(db, "SNP", df)
#    invisible(TRUE)
#}

library(dplyr)

snpCollapsedView <- function(snp) {
    stopifnot(length(snp) >= 1)
    stopifnot(is(snp, "GRanges"))
    stopifnot(all(width(snp) == 1))
    
    snpdf <- tibble::as_tibble(as.data.frame(snp))
    snpdf[["width"]] <- NULL
    
    grouped <- dplyr::group_by(snpdf, seqnames, start, end, strand, REF)
    ans <- dplyr::summarize(grouped, ALL = paste0(REF, ALT, collapse = ""))
    ans <- dplyr::ungroup(ans)
    ans$NS <- unique(snpdf$NS)
    ans$AN <- unique(snpdf$AN)
    
    # Summarize AC
    summarize_AC <- function(char) {
        tmp <- dplyr::select(snpdf, seqnames, start, end, strand, ALT, AC)
        tmp <- dplyr::filter(tmp, ALT == !!char)
        new_colname <- paste0("AC", "_", char)
        tmp <- dplyr::rename(tmp, !!new_colname := AC)
        tmp <- dplyr::select(tmp, -ALT)
        tmp
    }
    # Join the columns
    ## FIXME: I ignored the `Het` and `Hom` columns.
    ans <- local({
        nrow_before <- nrow(ans)
        ans <- dplyr::left_join(ans, summarize_AC("A"))
        ans <- dplyr::left_join(ans, summarize_AC("C"))
        ans <- dplyr::left_join(ans, summarize_AC("G"))
        ans <- dplyr::left_join(ans, summarize_AC("T"))
        stopifnot(nrow(ans) == nrow_before)
        ans
    })
    
    ans$ALL <- mergeIUPACLetters(ans$ALL)
    ans
}

suppressWarnings({
    chain.panTro5 <- import.chain(here("raw_data/chainfiles/hg38ToPanTro5.over.chain"))
    chain.gorGor5 <- import.chain(here("raw_data/chainfiles/hg38ToGorGor5.over.chain"))
    chain.ponAbe2 <- import.chain(here("raw_data/chainfiles/hg38ToPonAbe2.over.chain"))
    chain.rheMac8 <- import.chain(here("raw_data/chainfiles/hg38ToRheMac8.over.chain"))
})

lift_over <- function(snpdf, bsgenome, chain) {
    stopifnot(is.data.frame(snpdf))
    stopifnot(is(chain, "Chain"))
    
    gr <- as(snpdf, "GPos")
    mcols(gr) <- NULL
    
    grl.target <- rtracklayer::liftOver(gr, chain)
    gr$target.loc <- grl.target
    
    queryString <- function(grl, bsgenome) {
        gr <- unlist(unname(grl))
        dna <- DNAStringSet(Views(bsgenome, gr))
        relist(dna, unname(grl))
    }
    
    gr$target.DNA <- queryString(grl.target, bsgenome)
    
    # Unlist the DNAStringSetList
    unlist_dna <- function(dsslist) {
        len_ori <- length(dsslist)
        dsslist[lengths(dsslist) == 0] <- DNAStringSetList(DNAStringSet("-"))
        ans <- unlist(dsslist)
        stopifnot(length(ans) == len_ori)
        ans
    }
    gr$target.DNA <- unlist_dna(gr$target.DNA)
    
    ans <- mcols(gr)
    
    colnames(ans) <- c(paste0("loc_", providerVersion(bsgenome)),
                       paste0("dna_", providerVersion(bsgenome)))
    ans
}

lift_all <- function(snpdf) {
    d1 <- (lift_over(snpdf, BSgenome.Ptroglodytes.UCSC.panTro5, chain.panTro5))
    d2 <- (lift_over(snpdf, BSgenome.GGorilla.UCSC.gorGor5, chain.gorGor5))
    d3 <- (lift_over(snpdf, BSgenome.PAbelii.UCSC.ponAbe2, chain.ponAbe2))
    d4 <- (lift_over(snpdf, BSgenome.Mmulatta.UCSC.rheMac8, chain.rheMac8))
    # I dropped the location information
    ans <- cbind(d1, d2, d3, d4)[,c(2,4,6,8)]
    tibble::as_tibble(as.data.frame(ans))
}



## Phangorn staff
library(phangorn)
TREE = local({
    tree <- read.tree(
        text = "((((hg:6.4[&&NHX:dist=6.4:name=hg:support=1.0],panTro:6.4[&&NHX:dist=6.4:name=panTro:support=1.0])1:2.21[&&NHX:dist=2.21:name=:support=1.0],gorGor:8.61[&&NHX:dist=8.61:name=gorGor:support=1.0])1:6.59[&&NHX:dist=6.59:name=:support=1.0],ponAbe:15.2[&&NHX:dist=15.2:name=ponAbe:support=1.0])1:12.9[&&NHX:dist=12.9:name=:support=1.0],rheMac:28.1[&&NHX:dist=28.1:name=rheMac:support=1.0]);"
        #"((((hg:6.4[&&NHX:dist=6.4:name=hg:support=1.0],panTro:6.4[&&NHX:dist=6.4:name=panTro:support=1.0])1:2.21[&&NHX:dist=2.21:name=hg.panTro:support=1.0],gorGor:8.61[&&NHX:dist=8.61:name=gorGor:support=1.0])1:6.59[&&NHX:dist=6.59:name=hg.panTro.ponAbe:support=1.0],ponAbe:15.2[&&NHX:dist=15.2:name=ponAbe:support=1.0])1:12.9[&&NHX:dist=12.9:name=root:support=1.0],rheMac:28.1[&&NHX:dist=28.1:name=rheMac:support=1.0]);"
    )
    #tree$node.label <- NULL
    tree$node.label <- c("root", "ponAbe-else", "hg-pan-gor", "hg-panTro")
    tree
})
#getMRCA(TREE, c("hg", 'panTro'))
#allDescendants(TREE)
if (FALSE)
    plot(TREE, show.node.label = TRUE)

COST_MAT <- function(ratio) {
    levels <- c("a", "c", "g", "t")
    transversion_cost <- ratio
    cost <- matrix(transversion_cost, 4, 4, dimnames = list(levels, levels))
    cost <- cost - transversion_cost * diag(4)
    cost
    transitions <- rbind(
        c("a", "g"),
        c("g", "a"),
        c("c", "t"),
        c("t", "c")
    )
    cost[transitions] <- 1
    #unname(cost)
    cost
}
COST_MAT(1.2)

fit_anc <- function(phyd) {
    TREE
    stopifnot(identical(names(phyd), TREE$tip.label))
    cost <- COST_MAT(1.2)
    stopifnot(identical(rownames(cost), attr(phyd, "levels")))
    stopifnot(identical(colnames(cost), attr(phyd, "levels")))
    
    ans <- ancestral.pars(
        tree = TREE,
        data = phyd,
        return = "phyDat",
        cost = cost,
        type = "MPR"
    )
    ans <- as.character(ans)
    stopifnot(identical(rownames(ans)[seq(Ntip(TREE))], TREE$tip.label))
    rownames(ans)[seq(Ntip(TREE) + 1, Ntip(TREE) + Nnode(TREE))] <- TREE$node.label
    ans <- t(ans)
    ans <- as.data.frame(ans, stringsAsFactors = FALSE)
    ans
}


## ======= Starting the tasks ===========

slidingRanges <- slidingGenome()
#slidingRanges <- sample(slidingGenome(chr = "chr8"), 10)

# #### ## Chr8
#set.seed(5)
#    slidingRanges <- slidingRanges[seqnames(slidingRanges) == "chr8"]
#    slidingRanges <- slidingRanges[sort(sample(1:length(slidingRanges), 10))]
#    slidingRanges
# #### ##

#rg <- slidingRanges[1]

RES_DIR = "results"
if (!dir.exists(RES_DIR))
    dir.create(RES_DIR, recursive = TRUE)

LOG_FILE <- "snp_db.log"

main_run <- function(rg) {
    
    identity_name <- sprintf("%s_%s", as.character(rg),
                             format(Sys.time(), "%Y-%m-%d_%H:%M"))
    output_path <- file.path(RES_DIR, paste0(identity_name, ".tsv.gz"))
    
    cat(sprintf("> Running %s\n", as.character(rg)), file = LOG_FILE, append = TRUE)
    
    snp <- extractSNPFromVcf(vcf = VCF_FILE_LOC, rg)
    colsnp <- snpCollapsedView(snp)
    lift_res <- lift_all(colsnp)
    lift_res <- cbind(colsnp, lift_res)
    phyd <- dplyr::select(
        lift_res,
        hg = ALL,
        panTro = dna_panTro5,
        gorGor = dna_gorGor5,
        ponAbe = dna_ponAbe2,
        rheMac = dna_rheMac8
    )
    phyd <- phyDat(t(as.matrix(phyd)))
    stopifnot(identical(names(phyd), TREE$tip.label))
    
    anc_res <- fit_anc(phyd)
    anc_res
    
    res <- dplyr::select(lift_res, - starts_with("dna"))
    res <- cbind(res, anc_res)
    res
    ## res is a data frame that we want to store
    
    readr::write_tsv(res, output_path)
    cat(sprintf("< Finished %s\n", as.character(rg)), file = LOG_FILE, append = TRUE)
}

bplapply(seq_along(slidingRanges), BPPARAM = MulticoreParam(),
    function(i) {
        main_run(slidingRanges[i])
    }
)

