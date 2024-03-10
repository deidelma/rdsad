"""

rds_to_h5ad.py

Script to process .rds files containing raw 10x data to h5ad files.

Designed to handle .rds files that are not compatible with anndata.read_10x_mtx
because they do not represent Seurat object. 
Works by extracting the information from the object assuming that there is
valid matrix, barcode, and feature information.

Copyright (c) 2024 David Eidelman, MIT License
"""
import scipy.io
import pandas as pd
import rpy2.robjects as robjects
from anndata import AnnData
from pathlib import Path

RDS_FILE = Path("data/JK06.rds") # need to replace this with a command line argument
DATA_DIR = Path("data")  # need to make this modifieable

robjects.r("library(Seurat, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)")
robjects.r('library(SeuratData, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)')
robjects.r('library(SeuratDisk, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)')
robjects.r('library(fastMatMR, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)')

read_rds = robjects.r("readRDS")
rds_data = read_rds(RDS_FILE.as_posix())

seurat_object = robjects.r["CreateSeuratObject"](rds_data)
mtx = robjects.r["GetAssayData"](seurat_object, assay='RNA', layer='counts')
robjects.r["write_fmm"](mtx, (DATA_DIR / Path("matrix.mtx")).as_posix())

rownames = robjects.r["rownames"]
barcodes = rownames(seurat_object.slots["meta.data"])

x = tuple(seurat_object.slots["assays"].items())[0][1]
features = rownames(x.slots["features"])

write_table = robjects.r["write.table"]
write_table(barcodes, 
            (DATA_DIR / Path('barcodes.tsv')).as_posix(),
            quote = False, row_names=False, col_names=False)
write_table(features, 
            (DATA_DIR / Path('features.tsv')).as_posix(),
            quote = False, row_names=False, col_names=False)

mtx = scipy.io.mmread(DATA_DIR / Path("matrix.mtx"))
barcodes = pd.read_csv(DATA_DIR / Path("barcodes.tsv"), header=None)
features = pd.read_csv(DATA_DIR / Path("features.tsv"), header=None)

adata = AnnData(mtx.tocsr().T)
adata.obs = barcodes
adata.obs_names = barcodes[0]
adata.obs.drop(adata.obs[[0]], axis=1, inplace=True)
adata.var = features
adata.var_names = features[0]
adata.var.drop(adata.var[[0]], axis=1, inplace=True)

adata.write_h5ad(DATA_DIR / "jk06.h5ad")

