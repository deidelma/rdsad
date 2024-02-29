#!/usr/bin/env Rscript

#
# rdsad-convert.R
#
# this script converts .rds files into .h5ad files using the approach outlined
# by Z.Q. Fang (https://zqfang.github.io/2020-04-28-seurat2scanpy).
# However, instead of using the method based on MuDataSeurat, which doesn't seem to work
# as of 2024-02-29, it uses the SeuratDisk method recommended by tehe SatijaLab, who
# are the authors of the Seurat system.
#
# Copyright (c) 2024 by David Eidelman. All rights reserved.  MIT License.
#
#

# to call this script you need to pass in the filename of the file to be converted 
# as a parameter to the script.
#
# For example:
#
#  ./rdsad-convert.R myfile.rds
#
# the output will be a file called myfile.h5ad
#
library(Seurat)
library(SeuratDisk)
library(stringr)
#
# check command line
#
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

#
# convert factor to character
#
i <- sapply(data@meta.data, is.factor)
data@meta.data[i] <- lapply(data@meta.data[i], as.character)

#
# convert file
#
fname <- stringr::str_replace(args[1], ".rds", ".h5seurat")
message(sprintf("Creating intermediate file: %s", fname))
SaveH5Seurat(data, filename=fname, overwrite=TRUE)
message("Converting file...")
Convert(fname, dest="h5ad", assay="RNA", overwrite=TRUE)
message("File successfully converted.")


