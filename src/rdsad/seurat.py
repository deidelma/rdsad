"""

seurat.py

Module to handle manipulation of Seurat and related files.

Based on rpy2

"""

import tempfile
from os import PathLike
from pathlib import Path
from typing import Any, Callable
import rpy2.robjects as robjects
import scipy.io
from scipy.sparse import coo_matrix
from anndata import AnnData
import pandas as pd


class LoadingFailed(Exception):
    """
    Raised when unable to successfully load R or its libraries
    """

    ...


def load_library(library_name: str, loadr: Callable[[str], Any] = robjects.r) -> None:
    cmd = f"library({library_name}, "
    cmd += "quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)"
    return loadr(cmd)


def load_R_libraries() -> None:
    robjects.r("library(Seurat, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)")
    # robjects.r("library(Seurat, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)")
    return
    load_library("Seurat")
    load_library("SeuratData")
    load_library("SeuratDisk")
    load_library("fastMatMR")


def is_instance(obj: Any, type_str: str) -> bool:
    f = robjects.r("is")
    result = f(obj, type_str)  # type: ignore
    return result[0]


def convert_to_seurat(obj: Any) -> Any:
    # if obj is a list, extract the first element
    if robjects.r("class")(obj)[0] == "list":  # type: ignore
        obj = obj.rx2(1)
    # if obj is a matrix, convert it to a SeuratObject
    if is_instance(obj, "dgCMatrix"):
        # need to convert it to a seurat object
        return robjects.r("CreateSeuratObject")(obj)  # type: ignore
    elif is_instance(obj, "Seurat"):
        # if its a Seurat object, just return it
        return obj
    else:
        # not a valid object for translation to Seurat
        raise ValueError("Cannot be interpreted as a Seurat object.")


def read_rds_file(filename: PathLike | str):
    if isinstance(filename, Path):
        filename = filename.as_posix()
    read_rds: Callable[[str], Any] = robjects.r("readRDS")  # type: ignore
    assert isinstance(filename, str)
    return read_rds(filename)


def save_10x_data(seurat_object, data_dir: Path) -> Any:
    mtx = robjects.r["GetAssayData"](seurat_object, assay="RNA", layer="counts")  # type: ignore
    robjects.r("fastMatMR::write_fmm")(mtx, (data_dir / Path("matrix.mtx")).as_posix())  # type: ignore
    # robjects.r["fastMatMR::write_fmm"](mtx, (data_dir / Path("matrix.mtx")).as_posix())  # type: ignore

    rownames = robjects.r["rownames"]
    barcodes = rownames(seurat_object.slots["meta.data"])  # type: ignore

    x = tuple(seurat_object.slots["assays"].items())[0][1]
    features = rownames(x.slots["features"])  # type: ignore

    write_table = robjects.r["write.table"]
    write_table(
        barcodes,
        (data_dir / Path("barcodes.tsv")).as_posix(),
        quote=False,
        row_names=False,
        col_names=False,
    )  # type: ignore
    write_table(
        features,
        (data_dir / Path("features.tsv")).as_posix(),
        quote=False,
        row_names=False,
        col_names=False,
    )  # type:ignore


def read_10x_data(data_dir: Path) -> tuple[coo_matrix, pd.DataFrame, pd.DataFrame]:
    mtx = scipy.io.mmread(data_dir / Path("matrix.mtx"))
    barcodes = pd.read_csv(data_dir / Path("barcodes.tsv"), header=None, sep="\t")
    features = pd.read_csv(data_dir / Path("features.tsv"), header=None, sep="\t")
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


def seurat_to_h5ad(obj, output_file: PathLike) -> None:
    if not is_instance(obj, "Seurat"):
        raise ValueError("Attempt to convert non-Seurat object to h5ad format.")
    with tempfile.TemporaryDirectory() as tmp_dir:
        dir_path = Path(tmp_dir)
        save_10x_data(obj, dir_path)
        mtx, barcodes, features = read_10x_data(dir_path)
        adata = create_adata(mtx=mtx, barcodes=barcodes, features=features)
        print("adata created")
        adata.write_h5ad(output_file)
