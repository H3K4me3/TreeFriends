
library(here)

## Unzip the file
if (!file.exists(here("raw_data/GTEx_v7_Binaryfile_4919.txt"))) {
    system2("unzip", args = c(here("raw_data/GTEx_v7_Binaryfile_4919.txt.zip"),
                              "-d", here("raw_data")))
}
stopifnot(file.exists(here("raw_data/GTEx_v7_Binaryfile_4919.txt")))          

## Read data
gtex <- readr::read_tsv(here("raw_data/GTEx_v7_Binaryfile_4919.txt"))
gtex
