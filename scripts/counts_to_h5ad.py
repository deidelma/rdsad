"""
counts_to_h5ad.py

Takes matrix, barcodes, and gene_names to create an h5ad file.

Useful when scanpy.read_10x_mtx is not suitable.
"""
import pandas as pd
import scipy.io
from anndata import AnnData, read_h5ad
from pathlib import Path

mtx = scipy.io.mmread("matrix.mtx")
# mtx = mtx.tocsr().T  # type: ignore
mtx = mtx.tocsr().T  # type: ignore

barcodes = pd.read_csv("barcodes.tsv", header=None)
features = pd.read_csv("features.tsv", header=None)

adata = AnnData(mtx)  # type:ignore
print("created adata")
adata.var = features  # barcodes
adata.var_names = adata.var[0]  # type: ignore
adata.var.drop(adata.var[[0]], axis=1, inplace=True)  # type:ignore
print("barcodes")
adata.obs = barcodes  # features
adata.obs_names = adata.obs[0]  # type:ignore
adata.obs.drop(adata.obs[[0]], axis=1, inplace=True)  # type: ignore
print("genes")
adata.write_h5ad(Path("success.h5ad"))
print("written")
adata = read_h5ad("success.h5ad")
print("read")
print(f"{adata}")
print(f"{adata.X}")
print(f"adata.var: {adata.var}")
print(f"adata.obs: {adata.obs}")

# sc.read_10x_mtx()
