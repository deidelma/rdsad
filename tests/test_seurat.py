from pathlib import Path
import tempfile
import sys
import os
from typing import Any
from rdsad import seurat
import rpy2.robjects as robjects
import pytest
import scanpy as sc


@pytest.fixture
def test_rds_file():
    dir_path = Path("tests/data")
    if not dir_path.exists():
        Path.mkdir(dir_path)

    rds_path = dir_path / Path("x.rds")
    filename = rds_path.with_suffix(".rds").as_posix()
    if not rds_path.exists():
        x = robjects.r("c(1, 2, 3)")
        save_rds = robjects.r("saveRDS")
        save_rds(x, filename)
    yield filename


def test_hello():
    assert True


def test_load_library() -> None:
    def loadr(x: str) -> Any:
        "dummy loader to be injected for testing purposes only"
        return x

    library_name = "Seurat"
    value = seurat.load_library(library_name, loadr=loadr)
    assert isinstance(value, str)
    assert "Seurat" in value
    assert (
        value
        == f"library({library_name}, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)"
    )


def test_load_R_libraries() -> None:
    # first ensure executing in right environment
    conda_prefix = os.environ["CONDA_PREFIX"]
    assert conda_prefix in ["rdsad", "/home/meakins/miniconda3"]
    # now try loading the libraries
    try:
        seurat.load_R_libraries()
    except Exception as e:
        print(f"ERROR:{e}", file=sys.stderr)
        assert False
    assert True


def test_read_rds(test_rds_file) -> None:
    seurat.load_R_libraries()
    # test reading generally
    read_rds = robjects.r("readRDS")
    x = read_rds(test_rds_file)  # type: ignore
    assert x is not None

    # ensure that read_rds_file detects invalid files
    with pytest.raises(ValueError):
        obj = seurat.read_rds_file("tests/data/x.rds")
        seurat.convert_to_seurat(obj)


@pytest.mark.skipif(
    not Path("tests/data/JK06.rds").exists(),
    reason="Need to have test JK file in tests/data",
)
def test_read_rds_file() -> None:
    seurat.load_R_libraries()
    file_path = Path("tests/data/JK06.rds")
    obj = seurat.read_rds_file(file_path)
    assert obj is not None

    seurat_obj = seurat.convert_to_seurat(obj)
    assert seurat_obj is not None
    assert seurat.is_instance(seurat_obj, "Seurat")


@pytest.mark.skipif(
    not Path("tests/data/pbmc.rds").exists(),
    reason="Need to have pbmc.rds file in tests/data",
)
def test_read_canonical_seurat_file() -> None:
    seurat.load_R_libraries()
    file_path = Path("tests/data/pbmc.rds")
    obj = seurat.read_rds_file(file_path)
    assert obj is not None
    assert seurat.is_instance(obj, "Seurat")

@pytest.mark.skipif(
    not Path("tests/data/pbmc.rds").exists(),
    reason="Need to have pbmc.rds file in tests/data",
)
def test_save_10x_data() -> None:
    seurat.load_R_libraries()
    file_path = Path("tests/data/pbmc.rds")
    obj = seurat.read_rds_file(filename=file_path)
    assert seurat.is_instance(obj, "Seurat")
    with tempfile.TemporaryDirectory() as tmp_dir:
    # tmp_dir = Path("tests/data/crap")
    # tmp_dir.mkdir()
        dir_path = Path(tmp_dir)
        seurat.save_10x_data(obj, dir_path)
        matrix_path = Path(dir_path / "matrix.mtx")
        features_path = Path(dir_path / "features.tsv")
        barcodes_path = Path(dir_path / "barcodes.tsv")
        assert matrix_path.exists()
        assert barcodes_path.exists()
        assert features_path.exists()
    # tmp_dir.rmdir()

        mtx, barcodes, features = seurat.read_10x_data(dir_path)
        assert mtx is not None
        assert barcodes is not None
        assert features is not None

        adata = seurat.create_adata(mtx, barcodes, features)
        assert adata is not None

@pytest.mark.skipif(
    not Path("tests/data/JK06.rds").exists(),
    reason="Need to have JK06.rds file in tests/data",
)
def test_seurat_to_h5ad() -> None:
    seurat.load_R_libraries()
    file_path = Path("tests/data/JK06.rds")
    matrix_obj = seurat.read_rds_file(file_path)
    obj = seurat.convert_to_seurat(matrix_obj)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        h5ad_file_path = (temp_path / file_path.name).with_suffix(".h5ad")
        seurat.seurat_to_h5ad(obj, h5ad_file_path)
        assert h5ad_file_path.exists() and h5ad_file_path.is_file()

        adata = sc.read_h5ad(h5ad_file_path)
        assert adata is not None
    