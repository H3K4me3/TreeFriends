
# library(here)
# 
# ## Unzip the file
# if (!file.exists(here("raw_data/GTEx_v7_Binaryfile_4919.txt"))) {
#     system2("unzip", args = c(here("raw_data/GTEx_v7_Binaryfile_4919.txt.zip"),
#                               "-d", here("raw_data")))
# }
# stopifnot(file.exists(here("raw_data/GTEx_v7_Binaryfile_4919.txt")))          
# 
# ## Read data
# gtex <- readr::read_tsv(here("raw_data/GTEx_v7_Binaryfile_4919.txt"))
# gtex
#


library(here)
library(readr)
ancres_data <- readr::read_tsv(here("results/GTEx_v7_Binaryfile_4919_ancestral_allele.txt"))
ancres_data

gtex_data <- readr::read_tsv(here("raw_data/GTEx_v7_Binaryfile_4919.txt"))
gtex_data

nrow(gtex_data) == nrow(ancres_data)

mean(ancres_data$ancestry1_panTro, na.rm = TRUE)
mean(ancres_data$ancestry1_panTro[which(!is.na(gtex_data$pval_Brain_Cortex))], na.rm = TRUE)
mean(ancres_data$ancestry1_panTro[which(!is.na(gtex_data$pval_Vagina))], na.rm = TRUE)

library(ggplot2)
