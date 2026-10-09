"""Read the compact public records and freshly generated, uncompressed records."""

import builtins
import glob
import gzip
from pathlib import Path


def open_result(file, mode="r", *args, **kwargs):
    if "r" in mode and isinstance(file, (str, Path)):
        path = Path(file)
        if not path.exists() and Path(str(path) + ".gz").exists():
            path = Path(str(path) + ".gz")
        if path.suffix == ".gz":
            return gzip.open(path, mode if "b" in mode or "t" in mode else mode + "t", *args, **kwargs)
        file = path
    return builtins.open(file, mode, *args, **kwargs)


def glob_results(pattern):
    native = glob.glob(pattern)
    compressed = glob.glob(pattern + ".gz")
    return native + [p for p in compressed if p[:-3] not in native]


def result_exists(path):
    return Path(path).exists() or Path(str(path) + ".gz").exists()
