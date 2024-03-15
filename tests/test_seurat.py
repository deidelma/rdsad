from pathlib import Path
import sys, os
from typing import Any
from rdsad import seurat
import rpy2.robjects as robjects
import pytest
import sh


def test_hello():
    assert True


def test_load_library() -> None:
    def loadr(x: str) -> Any:
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
    assert "rdsad" in os.environ["CONDA_PREFIX"]
    # now try loading the libraries
    try:
        seurat.load_R_libraries()
    except Exception as e:
        print(f"ERROR:{e}", file=sys.stderr)
        assert False
    assert True


def test_read_rds() -> None:
    seurat.load_R_libraries()
    # test reading generally
    read_rds = robjects.r("readRDS")
    x = read_rds('tests/data/x.rds') # type: ignore
    assert x is not None

    # ensure that read_rds_file detects invalid files
    with pytest.raises(ValueError):
        obj = seurat.read_rds_file("tests/data/x.rds")
        seurat.convert_to_seurat(obj)

@pytest.mark.skipif(not Path("tests/data/JK06.rds").exists(),
                    reason="Need to have test JK file in tests/data")
def test_read_rds_file() -> None:
    seurat.load_R_libraries()
    file_path = Path("tests/data/JK06.rds")
    obj = seurat.read_rds_file(file_path)
    assert obj is not None

    seurat_obj = seurat.convert_to_seurat(obj)
    assert seurat_obj is not None
