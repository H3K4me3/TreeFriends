
library(here)
library(magrittr)
library(VariantAnnotation)
library(rtracklayer)

if (!requireNamespace('BSgenome.Ptroglodytes.UCSC.panTro5')) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Ptroglodytes.UCSC.panTro5")
}
library(BSgenome.Ptroglodytes.UCSC.panTro5)


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

topmed.chr21 <- readVcf("rawdata/chr21.TOPMed_freeze5_hg38_dbSNP.recode.vcf.gz")
gr.hg38 <- rowRanges(topmed.chr21)

#### Liftover to panTro5  ###-------------------------------------------------------------------

grl.PanTro5 <- local({
    mcols(gr.hg38) <- NULL
    liftOver(gr.hg38, chain.PanTro5)
})
grl.PanTro5

dna.PanTro5 <- DNAStringSet(Views(BSgenome.Ptroglodytes.UCSC.panTro5, unlist(grl.PanTro5)))
names(dna.PanTro5) <- NULL
dna.PanTro5 <- relist(dna.PanTro5, grl.PanTro5)

stopifnot(identical(names(dna.PanTro5), names(gr.hg38)))
stopifnot(identical(names(dna.PanTro5), names(topmed.chr21)))

gr.hg38$DNA.PanTro5 <- dna.PanTro5

gr.hg38

saveRDS(gr.hg38, file = file.path(here(), "tmp_results/chr21_hg38.rds"))







