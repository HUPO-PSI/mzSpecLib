import sqlite3
import zlib
from dataclasses import dataclass

from typing import Iterator, List, Mapping, Tuple, Iterable, Type

import numpy as np

from pyteomics import proforma

from mzlib import annotation
from mzlib.analyte import FIRST_ANALYTE_KEY, FIRST_INTERPRETATION_KEY, Analyte
from mzlib.spectrum import Spectrum, SPECTRUM_NAME, CHARGE_STATE
from mzlib.attributes import AttributeManager, Attributed

from mzlib.backends.base import SpectralLibraryBackendBase, FORMAT_VERSION_TERM, DEFAULT_VERSION

from mzlib.index.base import IndexBase


def _compress_array(array: np.ndarray, dtype: str) -> bytes:
    """Compress the array to the EncyclopeDIA format."""
    packed = struct.pack(">" + (dtype * len(array)), *array)
    compressed = zlib.compress(packed, 9)
    return compressed


def _extract_array(byte_array: bytes, type_str="d") -> np.ndarray:
    dtype = np.dtype(type_str)
    decompressed = zlib.decompress(byte_array, 32)
    decompressed_length = len(decompressed) // dtype.itemsize
    unpacked = struct.unpack(">" + (type_str * decompressed_length), decompressed)
    return np.array(unpacked, dtype=dtype)


@dataclass
class EncyclopediaIndexRecord:
    number: int
    precursor_mz: float
    precursor_charge: int
    peptide: str


