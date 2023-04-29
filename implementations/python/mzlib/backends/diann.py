import csv
import io

from typing import List, Tuple, Dict, Iterator, Any

from pyteomics import proforma, mass

from mzlib import annotation
from mzlib.backends.base import SpectralLibraryBackendBase
from mzlib.backends.utils import open_stream
from mzlib.spectrum import Spectrum, SPECTRUM_NAME


def rewrite_unimod_peptide_as_proforma(sequence: str) -> str:
    return sequence.replace("(", '[').replace(')', ']').replace("UniMod", "UNIMOD")


CHARGE_STATE = "MS:1000041|charge state"
SELECTED_ION_MZ = "MS:1003208|experimental precursor monoisotopic m/z"
SOURCE_FILE = "MS:1003203|constituent spectrum file"
STRIPPED_PEPTIDE_TERM = "MS:1000888|stripped peptide sequence"
PROFORMA_PEPTIDE_TERM = "MS:1003169|proforma peptidoform sequence"


class DiaNNTSVSpectralLibrary(SpectralLibraryBackendBase):
    format_name = "dia-nn.tsv"

    def __init__(self, filename: str, index_type=None, **kwargs):
        if index_type is None:
            index_type = self.has_index_preference(filename)
        super().__init__(filename)
        self.filename = filename
        self._headers = None
        self._read_header_line()
        self.index, was_initialized = index_type.from_filename(filename)
        if not was_initialized:
            self.create_index()

    def _spectrum_type(self):
        key = "MS:1003072|spectrum origin type"
        value = "MS:1003074|predicted spectrum"
        return key, value

    def _read_header_line(self):
        with open_stream(self.filename) as stream:
            reader = csv.reader(stream, delimiter='\t')
            headers = next(reader)
            stream.seek(0)
        self._headers = headers

    def create_index(self):
        with open(self.filename, 'rb') as stream:
            header = stream.readline()
            header_cols = header.split(b'\t')
            column_key = header_cols.index(b'transition_group_id')
            offset = stream.tell()

            line = stream.readline()
            tokens = line.split(b'\t')
            key = tokens[column_key]
            index = {key: offset}
            n = 0
            self.index.add(
                number=n,
                offset=offset,
                name=key.decode("utf8"),
                analyte=None
            )
            n += 1

            # To hold previous values
            last_offset = None
            last_key = None
            while line:
                tokens = line.split(b'\t')
                if tokens[column_key] != key:
                    key = tokens[column_key]
                    self.index.add(
                        number=n,
                        offset=offset,
                        name=key.decode("utf8"),
                        analyte=None
                    )
                    n += 1
                    last_offset = stream.tell()
                else:
                    offset = stream.tell()
                line = stream.readline()
        self.index.add(
            number=n,
            offset=last_offset,
            name=tokens[column_key].decode('utf8'),
            analyte=None
        )
        n += 1
        self.index.commit()
        return n

    def _open_reader(self, stream: io.TextIOBase):
        return csv.DictReader(stream, fieldnames=self._headers, delimiter='\t')

    def _parse_from_buffer(self, buffer: List[Dict[str, Any]], spectrum_index: int = None) -> Spectrum:
        spec = self._new_spectrum()
        descr = buffer[0]

        spec.add_attribute(SPECTRUM_NAME, descr['transition_group_id'])
        spec.add_attribute(SELECTED_ION_MZ, float(descr['PrecursorMz']))
        spec.add_attribute(CHARGE_STATE, int(descr['PrecursorCharge']))
        spec.add_attribute(SOURCE_FILE, descr['FileName'])
        spec.add_attribute(*self._spectrum_type())

        analyte = self._new_analyte('1')

        pf_seq = rewrite_unimod_peptide_as_proforma(descr['FullUniModPeptideName'])
        peptide = proforma.ProForma.parse(pf_seq)

        analyte.add_attribute(STRIPPED_PEPTIDE_TERM, descr['PeptideSequence'])
        analyte.add_attribute(PROFORMA_PEPTIDE_TERM, pf_seq)
        analyte.add_attribute("MS:1001117|theoretical mass", peptide.mass)
        analyte.add_attribute("MS:1000885|protein accession", descr['UniprotID'])

        spec.add_analyte(analyte)
        spec.peak_list = self._generate_peaks(buffer)
        spec.add_attribute("MS:1003059|number of peaks", len(spec.peak_list))

        if spectrum_index:
            spec.index = spectrum_index
        else:
            spec.index = -1
        spec.key = spec.index + 1
        return spec

    def _generate_peaks(self, batch: List[Dict[str, Any]]) -> List[Tuple[float, float, List[annotation.IonAnnotationBase], List]]:
        peaks = []
        for row in batch:
            mz = float(row['ProductMz'])
            intensity = float(row['LibraryIntensity'])

            series = row['FragmentType']
            ordinal = int(row['FragmentSeriesNumber'])
            charge = int(row['FragmentCharge'])

            annot = annotation.PeptideFragmentIonAnnotation(
                series, ordinal, charge=charge, mass_error=annotation.MassError(0, 'Da'))

            peak = [
                mz, intensity, [annot], []
            ]
            peaks.append(peak)
        return peaks

    def read(self):
        with open_stream(self.filename) as stream:
            stream.readline()
            reader = self._open_reader(stream)
            batch_iter = batch_rows(reader)
            for rows in batch_iter:
                spec = self._parse_from_buffer(rows)
                yield spec

    def _get_lines_for(self, offset: int) -> List[Dict[str, Any]]:
        with open_stream(self.filename, 'r') as infile:
            infile.seek(offset)
            reader = self._open_reader(infile)
            spectrum_buffer = next(batch_rows(reader))
            #### We will end up here if this is the last spectrum in the file
        return spectrum_buffer

    def get_spectrum(self, spectrum_number: int = None, spectrum_name: str = None) -> Spectrum:
        # keep the two branches separate for the possibility that this is not possible with all
        # index schemes.
        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError("Provide only one of spectrum_number or spectrum_name")
            offset = self.index.offset_for(spectrum_number)
        elif spectrum_name is not None:
            offset = self.index.offset_for(spectrum_name)
        buffer = self._get_lines_for(offset)
        spectrum = self._parse_from_buffer(buffer, spectrum_number)
        return spectrum


def batch_rows(iterator: Iterator[Dict[str, Any]]) -> Iterator[List[Dict[str, Any]]]:
    group_key = None
    group = []
    for row in iterator:
        key = row['transition_group_id']
        if group_key is None:
            group_key = key
            group.append(row)
        elif group_key == key:
            group.append(row)
        else:
            yield group
            group = [row]
            group_key = key
    if group:
        yield group
