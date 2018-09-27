
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
stopifnot(file.exists("raw_data/chainfiles/hg38ToGorGor5.over.chain"))
stopifnot(file.exists("raw_data/chainfiles/hg38ToPanTro5.over.chain"))
stopifnot(file.exists("raw_data/chainfiles/hg38ToPonAbe2.over.chain"))
stopifnot(tools::md5sum("raw_data/chainfiles/hg38ToGorGor5.over.chain") == "fa532f74ede70ccc724badcf3a65acaf")
stopifnot(tools::md5sum("raw_data/chainfiles/hg38ToPanTro5.over.chain") == "cf39846b96f245d1ff27942d1ab94461")
stopifnot(tools::md5sum("raw_data/chainfiles/hg38ToPonAbe2.over.chain") == "18089c8a1b07268f8547a9402ce4d3b1")








