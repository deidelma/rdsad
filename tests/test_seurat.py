import sys, os
from typing import Any
from rdsad import seurat

# def test_load_R_libraries():


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
    try:
        seurat.load_R_libraries()
    except Exception as e:
        print(f"ERROR:{e}", file=sys.stderr)
        assert False
    assert True
