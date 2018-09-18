
library(here)
library(magrittr)
library(VariantAnnotation)
library(rtracklayer)
library(readr)


if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38")) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Hsapiens.UCSC.hg38")
}
library(BSgenome.Hsapiens.UCSC.hg38)

if (!requireNamespace('BSgenome.Ptroglodytes.UCSC.panTro5')) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Ptroglodytes.UCSC.panTro5")
}
if (!requireNamespace('BSgenome.GGorilla.UCSC.gorGor5'))
    devtools::install_git('https://gitlab.com/Marlin-Na/BSgenome.GGorilla.UCSC.gorGor5.git')

library(BSgenome.Ptroglodytes.UCSC.panTro5)
#library(BSgenome.GGorilla.UCSC.gorGor5)


#### Download chain files   ###----------------------------------------------------------------
chainfiledir <- file.path(here(), 'chainfiles')

if (!dir.exists(chainfiledir))
    dir.create(chainfiledir)

# Chimpanzee
if (!file.exists(file.path(chainfiledir, 'hg38ToPanTro5.over.chain'))) {
    download.file('http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToPanTro5.over.chain.gz',
                  destfile = file.path(chainfiledir, 'hg38ToPanTro5.over.chain.gz'))
    R.utils::gunzip(file.path(chainfiledir, "hg38ToPanTro5.over.chain.gz"), remove = FALSE)
}

# Gorilla
if (!file.exists(file.path(chainfiledir, 'hg38ToGorGor5.over.chain'))) {
    download.file('http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToGorGor5.over.chain.gz',
                  destfile = file.path(chainfiledir, 'hg38ToGorGor5.over.chain.gz'))
    R.utils::gunzip(file.path(chainfiledir, "hg38ToGorGor5.over.chain.gz"), remove = FALSE)
}

chain.PanTro5 <- import.chain(file.path(chainfiledir, "hg38ToPanTro5.over.chain"))
chain.GorGor5 <- import.chain(file.path(chainfiledir, "hg38ToGorGor5.over.chain"))

chain.PanTro5
chain.GorGor5

# LiftOver



#### Read Vcf files   ###----------------------------------------------------------------------

#topmed.chr21 <- readVcf("rawdata/chr21.TOPMed_freeze5_hg38_dbSNP.recode.vcf.gz")
df.topmed.sample <- readr::read_tsv(file = "raw_data/sample_lines.tsv", col_names = FALSE)
colnames(df.topmed.sample) <- c("seqnames", "pos", "id", "REF", "ALT", "QUAL", "FILTER", "INFO",
                                "FORMAT", "TOPMED")
library(dplyr)
df.topmed.sample <- dplyr::select(df.topmed.sample, seqnames, pos, id, REF, ALT)
df.topmed.sample <- dplyr::rename(df.topmed.sample, start = pos)
df.topmed.sample <- dplyr::mutate(df.topmed.sample, width = nchar(REF))
df.topmed.sample$strand <- "*"
df.topmed.sample <- dplyr::mutate(df.topmed.sample, end = start + width - 1)
gr.topmed.sample <- as(df.topmed.sample, "GRanges")
gr.topmed.sample



#### Liftover to panTro5  ###-------------------------------------------------------------------

grl.PanTro5 <- local({
    mcols(gr.topmed.sample) <- NULL
    liftOver(gr.topmed.sample, chain.PanTro5)
})
grl.PanTro5
gr.topmed.sample$pos.PanTro5 <- grl.PanTro5

dna.PanTro5 <- DNAStringSet(Views(BSgenome.Ptroglodytes.UCSC.panTro5, unlist(grl.PanTro5)))
names(dna.PanTro5) <- NULL
dna.PanTro5 <- relist(dna.PanTro5, grl.PanTro5)

stopifnot(identical(names(dna.PanTro5), names(gr.topmed.sample)))

gr.topmed.sample$DNA.PanTro5 <- dna.PanTro5

gr.topmed.sample




gr.topmed.sample$DNA.PanTro5.expand <- ({
    grl <- gr.topmed.sample$pos.PanTro5
    start(grl) <- start(grl) - 3L
    end(grl) <- end(grl) + 3L
    
    dna.PanTro5.expand <- DNAStringSet(Views(BSgenome.Ptroglodytes.UCSC.panTro5, unlist(grl)))
    names(dna.PanTro5.expand) <- NULL
    dna.PanTro5.expand <- relist(dna.PanTro5.expand, grl)
    
    dna.PanTro5.expand
})

gr.topmed.sample$REF.expand <- {
    gr <- (gr.topmed.sample)
    mcols(gr) <- NULL
    start(gr) <- start(gr) - 3L
    end(gr) <- end(gr) + 3L
    
    ref.expand <- DNAStringSet(Views(BSgenome.Hsapiens.UCSC.hg38, gr))
    names(ref.expand) <- NULL
    ref.expand
}

fliping <- mcols(gr.topmed.sample)[c("REF.expand", "DNA.PanTro5.expand")]
fliping <- fliping[lengths(fliping$DNA.PanTro5.expand) == 1, ]
fliping$DNA.PanTro5.expand <- unlist(fliping$DNA.PanTro5.expand)
fliping$DNA.PanTro5.expand.reverse <- Biostrings::reverseComplement(fliping$DNA.PanTro5.expand)
fliping
fliping[,c(1,2)]

sum((fliping$REF.expand == fliping$DNA.PanTro5.expand) | (fliping$REF.expand == fliping$DNA.PanTro5.expand.reverse))
fliping[!((fliping$REF.expand == fliping$DNA.PanTro5.expand) | (fliping$REF.expand == fliping$DNA.PanTro5.expand.reverse)), ]



