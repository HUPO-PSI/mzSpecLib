import sqlite3
import zlib
from dataclasses import dataclass

from typing import Any, Iterator, List, Mapping, Tuple, Iterable, Type

import numpy as np

from pyteomics import proforma

from mzlib import annotation
from mzlib.analyte import FIRST_ANALYTE_KEY, FIRST_INTERPRETATION_KEY, Analyte, ProteinDescription
from mzlib.spectrum import Spectrum, SPECTRUM_NAME, CHARGE_STATE
from mzlib.attributes import AttributeManager, Attributed, Attribute

from mzlib.backends.base import SpectralLibraryBackendBase, FORMAT_VERSION_TERM, DEFAULT_VERSION

from mzlib.index.base import IndexBase


DECOY_SPECTRUM = "MS:1003192|decoy spectrum"
DECOY_PEPTIDE_SPECTRUM = "MS:1003195|unnatural peptidoform decoy spectrum"


def _decode_peaks(record: sqlite3.Row):
    raw_data = zlib.decompress(record['MassArray'])
    mass_array = np.frombuffer(raw_data, dtype='>d')
    raw_data = zlib.decompress(record['IntensityArray'])
    intensity_array = np.frombuffer(raw_data, dtype='>f')
    return mass_array, intensity_array


@dataclass
class EncyclopediaIndexRecord:
    number: int
    precursor_mz: float
    precursor_charge: int
    peptide: str


class EncyclopediaIndex(IndexBase):
    connection: sqlite3.Connection

    def __init__(self, connection):
        self.connection = connection

    def __getitem__(self, i):
        if isinstance(i, int):
            return self.search(i + 1)
        elif isinstance(i, slice):
            return [self.search(j + 1) for j in range(i.start or 0, i.stop or len(self), i.step or 1)]
        else:
            raise TypeError(f"Cannot index {self.__class__.__name__} with {i}")

    def _record_from(self, row: Mapping) -> EncyclopediaIndexRecord:
        peptide_sequence = row['PeptideModSeq']
        return EncyclopediaIndexRecord(row['rowid'], row['PrecursorMz'], row['PrecursorCharge'], peptide_sequence)

    def search(self, i):
        if isinstance(i, int):
            info = self.connection.execute("SELECT rowid, PrecursorMz, PrecursorCharge, PeptideModSeq FROM entries WHERE rowid = ?", (i, )).fetchone()
            return self._record_from(info)
        elif isinstance(i, str):
            raise NotImplementedError()

    def __iter__(self):
        return map(self._record_from, self.connection.execute("SELECT rowid, PrecursorMz, PrecursorCharge, PeptideModSeq FROM entries ORDER BY rowid").fetchall())

    def __len__(self):
        return self.connection.execute("SELECT count(rowid) FROM entries;").fetchone()[0]


class EncyclopediaSpectralLibrary(SpectralLibraryBackendBase):
    """Read EncyclopeDIA SQLite3 spectral library files."""

    connection: sqlite3.Connection

    file_format = "dlib"
    format_name = "encyclopedia"

    @classmethod
    def has_index_preference(cls, filename) -> Type[IndexBase]:
        return EncyclopediaIndex

    def __init__(self, filename: str, **kwargs):
        super().__init__(filename)
        self.connection = sqlite3.connect(filename)
        self.connection.row_factory = sqlite3.Row
        self.index = EncyclopediaIndex(self.connection)
        self.read_header()

    def read_header(self) -> bool:
        attribs = AttributeManager()
        attribs.add_attribute(FORMAT_VERSION_TERM, DEFAULT_VERSION)
        attribs.add_attribute("MS:1003207|library creation software", "EncyclopeDIA")
        self.attributes = attribs
        return True

    def _populate_analyte(self, analyte: Analyte, row: Mapping[str, Any]):
        """
        Fill an analyte with details describing a peptide sequence and inferring
        from context its traits based upon the assumptions EncyclopeDIA makes.

        EncyclopeDIA only stores modifications as delta masses.
        """
        peptide = proforma.ProForma.parse(row['PeptideModSeq'])
        analyte.add_attribute("MS:1003169|proforma peptidoform sequence", str(peptide))
        analyte.add_attribute("MS:1001117|theoretical mass", peptide.mass)
        analyte.add_attribute("MS:1000888|stripped peptide sequence", row['PeptideSeq'])
        analyte.add_attribute(CHARGE_STATE, row['PrecursorCharge'])

        cursor = self.connection.execute(
            "SELECT ProteinAccession, isDecoy FROM peptidetoprotein WHERE PeptideSeq = ?;", (row['PeptideSeq'], ))

        had_decoy = False
        for protrow in cursor:
            accession = protrow['ProteinAccession']
            is_decoy = bool(int(protrow['isDecoy']))
            had_decoy = had_decoy or is_decoy
            analyte.add_attribute_group([
                Attribute('MS:1000885|protein accession', accession)
            ])
        return had_decoy

    def get_spectrum(self, spectrum_number: int = None, spectrum_name: str = None):
        """
        Read a spectrum from the spectrum library.

        EncyclopeDIA does not support alternative labeling of spectra with a
        plain text name so looking up by `spectrum_name` is not supported.
        """
        if spectrum_number is None:
            raise ValueError("Only spectrum number queries are supported. spectrum_number must have an integer value")
        try:
            spectrum_number = int(spectrum_number)
        except (TypeError, ValueError):
            raise ValueError(f"spectrum_number must have an integer value, received {spectrum_number!r}") from None

        info = self.connection.execute("SELECT rowid, * FROM entries WHERE rowid = ?;", (spectrum_number, )).fetchone()
        spectrum = self._new_spectrum()
        spectrum.key = info['rowid']
        spectrum.index = info['rowid'] - 1
        spectrum.precursor_mz = info['PrecursorMz']
        try:
            spectrum.add_attribute("MS:1000894|retention time", info['RTInSeconds'] / 60.0)
        except KeyError:
            pass

        try:
            spectrum.add_attribute(
                "MS:1003203|constituent spectrum file",
                f"file://{info['SourceFile']}"
            )
        except KeyError:
            pass


        analyte = self._new_analyte(1)
        had_decoy = self._populate_analyte(analyte, info)
        if had_decoy:
            spectrum.add_attribute(DECOY_SPECTRUM, DECOY_PEPTIDE_SPECTRUM)

        spectrum.add_analyte(analyte)

        interp = self._new_interpretation(1)
        interp.add_analyte(analyte)
        spectrum.add_interpretation(interp)

        mz_array, intensity_array = _decode_peaks(info)
        n_peaks = len(mz_array)
        spectrum.add_attribute("MS:1003059|number of peaks", n_peaks)

        peak_list = []
        # EncyclopeDIA does not encode product ion identities
        for i, mz in enumerate(mz_array):
            row = (mz, intensity_array[i], [], [])
            peak_list.append(row)
        spectrum.peak_list = peak_list
        return spectrum

    def read(self) -> Iterator[Spectrum]:
        for rec in self.index:
            yield self.get_spectrum(rec.number)
