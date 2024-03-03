#!/usr/bin/env Rscript
#
# extract_items.R
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

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	stop("No filename provided.", call.=FALSE)
}

#
# load the file
#
message(sprintf("Processing %s", args[1]))
data <- readRDS(args[1])
message("File successfully read.")
