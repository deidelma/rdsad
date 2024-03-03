#!/usr/bin/env Rscript

#
# rdsad-install.R 
#
# This script loads the packages SeuratData and SeuratDisk from their github repositories,
# thus ensuring compatability with the current version of Seurat.
#
# Copyright (c) 2024, David Eidelman, MIT license.
#
library(Seurat)
print("Seurat successfully loaded")
# install using the default 
devtools::install_github("satijalab/seurat-data", upgrade=TRUE)
print("SeuratData successfully installed")
devtools::install_github("mojaveazure/seurat-disk", upgrade=TRUE)
print("SeuratDisk successfully installed")
