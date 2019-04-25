

suppressPackageStartupMessages({
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
    library(phangorn)
})

#options(verbose = TRUE)
options(mc.cores = 4)

source(here("lib/liftover.R"))


VCF_FILE_LOC = here("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz")
#DB_LOC = here("results/snpdb.sqlite")

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


suppressWarnings({
    chain.panTro5 <- import.chain(here("raw_data/chainfiles/hg38ToPanTro5.over.chain"))
    chain.gorGor5 <- import.chain(here("raw_data/chainfiles/hg38ToGorGor5.over.chain"))
    chain.ponAbe2 <- import.chain(here("raw_data/chainfiles/hg38ToPonAbe2.over.chain"))
    chain.rheMac8 <- import.chain(here("raw_data/chainfiles/hg38ToRheMac8.over.chain"))
})

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
#set.seed(4)
#slidingRanges <- sample(slidingGenome(chr = "chr8"), 10)

# #### ## Chr8
#set.seed(5)
#    slidingRanges <- slidingRanges[seqnames(slidingRanges) == "chr8"]
#    slidingRanges <- slidingRanges[sort(sample(1:length(slidingRanges), 10))]
#    slidingRanges
# #### ##

#rg <- slidingRanges[1]

write_log <- function(file, content) {
    date <- format(Sys.time(), "%Y-%m-%d_%H:%M")
    date <- paste0("@", date, " > ")
    cat(paste0(date, content, "\n"), file = file, append = TRUE)
}

RES_DIR = here("results/snptable")
if (!dir.exists(RES_DIR))
    dir.create(RES_DIR, recursive = TRUE)

main_run <- function(rg) {
    
    identity_name <- sprintf("%s_%s", as.character(rg),
                             format(Sys.time(), "%Y-%m-%d_%H:%M"))
    output_path <- file.path(RES_DIR, paste0(identity_name, ".tsv.gz"))
    LOG_FILE <- file.path(RES_DIR, paste0(identity_name, ".log"))
    
    write_log(LOG_FILE, "Start Running")
    write_log(LOG_FILE, "Running ReadVcf")
    
    snp <- extractSNPFromVcf(vcf = VCF_FILE_LOC, rg)
    
    if (length(snp) == 0) {
        write_log(LOG_FILE, "Skipping liftover etc..")
        
        write_log(LOG_FILE, "Running WriteData")
        
        # Write an empty data frame
        readr::write_tsv(tibble::tibble(), output_path)
        
        write_log(LOG_FILE, "Finished Running All")
        
        return()
    }
    
    write_log(LOG_FILE, "Running CollapsedView")
    
    colsnp <- snpCollapsedView(snp)
    
    write_log(LOG_FILE, "Running LiftAll")
    
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
    
    write_log(LOG_FILE, "Running FitAnc")
    
    anc_res <- fit_anc(phyd)
    anc_res
    
    res <- dplyr::select(lift_res, - dplyr::starts_with("dna"))
    res <- cbind(res, anc_res)
    res
    ## res is a data frame that we want to store
    
    write_log(LOG_FILE, "Running WriteData")
    
    readr::write_tsv(res, output_path)
    
    write_log(LOG_FILE, "Finished Running All")
    NULL
}

bplapply(seq_along(slidingRanges), BPPARAM = MulticoreParam(),
    function(i) {
        main_run(slidingRanges[i])
    }
)

