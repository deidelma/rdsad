"""

seurat.py

Module to handle manipulation of Seurat and related files.

Based on rpy2

"""
from os import PathLike
from pathlib import Path
from enum import Enum
from typing import Any, Callable
import rpy2.robjects as robjects


class LoadingFailed(Exception):
    """
    Raised when unable to successfully load R or its libraries
    """
    ...

def load_library(library_name: str, loadr: Callable[[str],Any] = robjects.r) -> None:
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


# def load_R_libraries() -> None:
#     """
#     Loads R and required libraries.

#     Raises LoadingFailed if an error is encountered.
#     """
#     try:
#         robjects.r(
#                 "library(Seurat, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)")
#         robjects.r(
#             "library(SeuratData, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)"
#         )
#         robjects.r(
#             "library(SeuratDisk, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)"
#         )
#         robjects.r("library(fastMatMR, quietly=TRUE, verbose= FALSE, warn.conflicts=FALSE)")
#     except Exception:
#         raise LoadingFailed("Unable to load R or its libraries")

def convert_to_seurat(obj: Any) -> Any:
    # if obj is a matrix, convert it to a SeuratObject
    get_desc: Callable[[Any],str] = robjects.r("str") # type: ignore
    desc = get_desc(obj)
    if "Formal class 'dgCMatrix'" in desc:
        # need to convert it to a seurat object
        seurat_obj = robjects.r("CreateSeuratObject")(obj) # type: ignore
    elif "Formal class 'Seurat'" in desc:
        seurat_obj = obj
    else:
        # not a valid object for translation to Seurat
        raise ValueError("'%s' cannot be interpreted as a Seurat object.", str(obj))
    
def read_rds_file(filename: PathLike| str):
    if isinstance(filename, Path):
        filename = filename.as_posix()
    read_rds: Callable[[str],Any] = robjects.r("readRDS") # type: ignore
    assert isinstance(filename,str)
    return read_rds(filename)