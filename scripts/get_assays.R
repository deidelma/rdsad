#!/usr/bin/env Rscript
#
# get_assays.R
#
# script to extract items from an .rds file containing a list of individual
# single cell RNA seq experiments and store them as .rds files
#
# based on https://stackoverflow.com/questions/59838033/saving-dataframes-stored-in-a-list-to-individual-files-in-r
#
# Copyright (c) 2024, David Eidelman. MIT License.
library(Seurat)
library(SeuratDisk)
library(stringr)
library(fastMatMR)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No filename provided.", call. = FALSE)
}

#
# load the file
#
message(sprintf("Processing %s", args[1]))
data <- readRDS(args[1])
message("File successfully read.")

if (is.list(data)) {
  data <- data[1]
}

data_seurat <- CreateSeuratObject(data)
mtx <- GetAssayData(data_seurat, assay = "RNA", layer = "counts")
write_fmm(mtx, "matrix.mtx")

barcodes <- rownames(data_seurat@meta.data)
write.table(barcodes,
  file = "barcodes.tsv",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

features <- rownames(data_seurat[["RNA"]]@features)
write.table(features,
  file = "features.tsv",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
