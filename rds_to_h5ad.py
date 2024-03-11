"""

rds_to_h5ad.py

Script to process .rds files containing raw 10x data to h5ad files.

Designed to handle .rds files that are not compatible with anndata.read_10x_mtx
because they do not represent Seurat object. 
Works by extracting the information from the object assuming that there is
valid matrix, barcode, and feature information.

Copyright (c) 2024 David Eidelman, MIT License
"""

import argparse
import sys
from typing import Any
import scipy.io
from scipy.sparse import coo_matrix
import pandas as pd
import rpy2.robjects as robjects
from anndata import AnnData
from pathlib import Path


def load_R_libraries() -> None:
    robjects.r("library(Seurat, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)")
    robjects.r(
        "library(SeuratData, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)"
    )
    robjects.r(
        "library(SeuratDisk, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)"
    )
    robjects.r("library(fastMatMR, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)")


def read_rds_file(filename: str) -> Any:
    read_rds = robjects.r("readRDS")
    rds_data = read_rds(filename)

    seurat_object = robjects.r["CreateSeuratObject"](rds_data)
    return seurat_object


def save_10x_data(seurat_object, data_dir: Path) -> Any:
    mtx = robjects.r["GetAssayData"](seurat_object, assay="RNA", layer="counts")
    robjects.r["write_fmm"](mtx, (data_dir / Path("matrix.mtx")).as_posix())

    rownames = robjects.r["rownames"]
    barcodes = rownames(seurat_object.slots["meta.data"])

    x = tuple(seurat_object.slots["assays"].items())[0][1]
    features = rownames(x.slots["features"])

    write_table = robjects.r["write.table"]
    write_table(
        barcodes,
        (data_dir / Path("barcodes.tsv")).as_posix(),
        quote=False,
        row_names=False,
        col_names=False,
    )
    write_table(
        features,
        (data_dir / Path("features.tsv")).as_posix(),
        quote=False,
        row_names=False,
        col_names=False,
    )


def read_10x_data(data_dir: Path) -> tuple[coo_matrix, pd.DataFrame, pd.DataFrame]:
    mtx = scipy.io.mmread(data_dir / Path("matrix.mtx"))
    barcodes = pd.read_csv(data_dir / Path("barcodes.tsv"), header=None)
    features = pd.read_csv(data_dir / Path("features.tsv"), header=None)
    return coo_matrix(mtx), barcodes, features


def create_adata(
    mtx: coo_matrix, barcodes: pd.DataFrame, features: pd.DataFrame
) -> AnnData:
    adata = AnnData(mtx.tocsr().T)
    adata.obs = barcodes
    adata.obs_names = barcodes[0]  # type:ignore
    adata.obs.drop(adata.obs[[0]], axis=1, inplace=True)  # type: ignore
    adata.var = features
    adata.var_names = features[0]  # type:ignore
    adata.var.drop(adata.var[[0]], axis=1, inplace=True)  # type: ignore
    return adata


def main() -> int:
    parser = argparse.ArgumentParser(description="extract 10x data from .rds files")
    parser.add_argument("filename")
    parser.add_argument("--datadir", default="data")
    args = parser.parse_args()
    filename = args.filename
    if not Path(filename).exists():
        print(f"Unable to find file: {filename}. Exiting.", file=sys.stderr)
        return 1
    data_dir = Path(args.datadir)
    if not data_dir.exists() or not data_dir.is_dir():
        print(f"{data_dir} is not a valid directory.  Exiting.")
        return 1
    print(f"Processing file: {args.filename}")
    load_R_libraries()
    print("R libraries loaded.")
    seurat_object = read_rds_file(filename=filename)
    print("Seurat object created")
    save_10x_data(seurat_object=seurat_object, data_dir=data_dir)
    print("10X data saved")
    mtx, barcodes, features = read_10x_data(data_dir)
    adata = create_adata(mtx=mtx, barcodes=barcodes, features=features)
    print("adata created")
    adata.write_h5ad(data_dir / "jk06.h5ad")
    print("h5ad file written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
