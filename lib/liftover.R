
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

snpCollapsedView <- function(snp) {
    stopifnot(length(snp) >= 1)
    stopifnot(is(snp, "GRanges"))
    stopifnot(all(width(snp) == 1))
    
    snpdf <- tibble::as_tibble(as.data.frame(snp))
    snpdf[["width"]] <- NULL
    
    grouped <- dplyr::group_by(snpdf, seqnames, start, end, strand, REF)
    ans <- dplyr::summarize(grouped, ALL = paste0(REF, ALT, collapse = ""))
    ans <- dplyr::ungroup(ans)
    #if (nrow(ans) == 0) {
    #    ans$ALL <- as.character(ans$ALL) # When snp has zero row, "ALL" will become a logical vector
    #}
    stopifnot(is.character(ans$ALL))
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
    
    #if (length(ans$ALL))
    #    ans$ALL <- mergeIUPACLetters(ans$ALL)
    ans$ALL <- mergeIUPACLetters(ans$ALL)
    ans
}

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
