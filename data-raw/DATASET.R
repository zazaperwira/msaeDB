## code to prepare `DATASET` dataset goes here

setwd("D:/1. SKRIPSI")
datamsaeDB <- read.csv("datamsaeDB2.csv",sep=";")


usethis::use_data(datamsaeDB, overwrite = TRUE)
