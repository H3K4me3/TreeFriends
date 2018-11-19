
library(here)

setwd(here())

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

### dbSNP package from Bioconductor  #####-------------------------------

if (!requireNamespace('SNPlocs.Hsapiens.dbSNP151.GRCh38'))
    BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")

### BSgenome packages from Bioconductor #####----------------------------

if (!requireNamespace('BSgenome.Hsapiens.UCSC.hg38'))
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

if (!requireNamespace('BSgenome.Ptroglodytes.UCSC.panTro5'))
    BiocManager::install("BSgenome.Ptroglodytes.UCSC.panTro5")

if (!requireNamespace('BSgenome.Mmulatta.UCSC.rheMac8'))
    BiocManager::install('BSgenome.Mmulatta.UCSC.rheMac8')

### Install self-made BSgenome packages #####---------------------------- 

if (!requireNamespace('BSgenome.GGorilla.UCSC.gorGor5'))
    devtools::install_url("https://github.com/H3K4me3/BSgenome.GGorilla.UCSC.gorGor5/releases/download/virtual-0.0.1/BSgenome.GGorilla.UCSC.gorGor5_0.0.1.tar.gz")

if (!requireNamespace('BSgenome.PAbelii.UCSC.ponAbe2'))
    devtools::install_url("https://github.com/H3K4me3/BSgenome.PAbelii.UCSC.ponAbe2/releases/download/virtual-0.0.1/BSgenome.PAbelii.UCSC.ponAbe2_0.0.1.tar.gz")

### Check raw data  #####------------------------------------------------

# File integrity
stopifnot(file.exists("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz"))
stopifnot(tools::md5sum("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz") == "773e9e97759a4a5b4555c5d7e1e14313")

# Create index
if (!file.exists("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz.tbi"))
    Rsamtools::indexTabix("raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz", format = "vcf")

# Check chain files
# The chain files are downloaded from ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
prepare_chainfile <- function(url, md5 = NULL) {
    gzpath <- here::here("raw_data/chainfiles", basename(url))
    path <- here::here("raw_data/chainfiles", sub("\\.gz$", "", basename(url)))
    if (!file.exists(path)) {
        download.file(url, gzpath)
        system2("gunzip", gzpath, wait = TRUE)
    }
    if (!is.null(md5))
        stopifnot(tools::md5sum(path) == md5)
    invisible(path)
}

prepare_chainfile(
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToGorGor5.over.chain.gz",
    "fa532f74ede70ccc724badcf3a65acaf"
)
prepare_chainfile(
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToPanTro5.over.chain.gz",
    "cf39846b96f245d1ff27942d1ab94461"
)
prepare_chainfile(
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToPonAbe2.over.chain.gz",
    "18089c8a1b07268f8547a9402ce4d3b1"
)
prepare_chainfile(
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToRheMac8.over.chain.gz",
    "da89c3353a70db359210ff7d45febf8d"
)






