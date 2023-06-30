import json
import os

from typing import List, Tuple, Dict, Iterator, Any, Deque, Union

from pyteomics import proforma

from mzlib import annotation
from mzlib.analyte import Analyte
from mzlib.backends.base import LIBRARY_NAME_TERM, _CSVSpectralLibraryBackendBase, FORMAT_VERSION_TERM, DEFAULT_VERSION
from mzlib.backends.utils import open_stream
from mzlib.spectrum import Spectrum, SPECTRUM_NAME


CHARGE_STATE = "MS:1000041|charge state"
SELECTED_ION_MZ = "MS:1003208|experimental precursor monoisotopic m/z"
SOURCE_FILE = "MS:1003203|constituent spectrum file"
STRIPPED_PEPTIDE_TERM = "MS:1000888|stripped peptide sequence"
PROFORMA_PEPTIDE_TERM = "MS:1003169|proforma peptidoform sequence"

CUSTOM_ATTRIBUTE_NAME = "MS:1003275|other attribute name"
CUSTOM_ATTRIBUTE_VALUE = "MS:1003276|other attribute value"

ID_SEP = "/"
_ID_SEP = ID_SEP.encode("ascii")

NO_LOSS = 'noloss'


def _rewrite_modified_peptide_as_proforma(sequence: str) -> str:
    if sequence.startswith("_"):
        sequence = sequence.strip("_")
    buffer: Deque[str] = Deque()
    last_paren = None
    for i, c in enumerate(sequence):
        if c == ']':
            if last_paren is not None:
                k = i - last_paren
                for j in range(k + 1):
                    buffer.pop()
                last_paren = None
                buffer.append(c)
        elif c == '(':
            last_paren = i
            buffer.append(c)
        else:
            buffer.append(c)
    return ''.join(buffer)


def _parse_value(value: str) -> Union[float, int, str, bool]:
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        lower = value.lower()
        if lower == "true":
            return True
        elif lower == "false":
            return False
        return value


class SpectronautTSVSpectralLibrary(_CSVSpectralLibraryBackendBase):
    format_name = "spectronaut.tsv"

    _custom_spectrum_keys = [
        "ExcludeFromAssay",
        "BGSInferenceId",
        "AllowForNormalization",
        "Workflow",
    ]

    _custom_analyte_keys = [
        "IsProteotypic",
        "FASTAName",
        "Database",
        "ProteinGroups",
    ]

    _key_columns = ['ModifiedPeptide', 'PrecursorCharge']

    _required_columns = ['PrecursorMz', 'PrecursorCharge',
                         'ModifiedPeptide', 'StrippedPeptide',
                         'FragmentMz', 'RelativeIntensity',
                         'FragmentType', 'FragmentNumber',
                         'FragmentCharge', 'FragmentLossType']

    def __init__(self, filename: str, index_type=None, **kwargs):
        super().__init__(filename, index_type=index_type, delimiter='\t', **kwargs)

    def _spectrum_type(self):
        key = "MS:1003072|spectrum origin type"
        value = "MS:1003073|observed spectrum"
        return key, value

    def read_header(self) -> bool:
        result = super().read_header()
        self.add_attribute(FORMAT_VERSION_TERM, DEFAULT_VERSION)
        if hasattr(self.filename, 'name'):
            name = self.filename.name.replace(".gz", '').rsplit('.', 1)[0].split(os.sep)[-1]
        else:
            name = self.filename.replace(".gz", '').rsplit(".", 1)[0].split(os.sep)[-1]
        self.add_attribute(LIBRARY_NAME_TERM, name)
        self.add_attribute("MS:1003207|library creation software", "MS:1001327|Spectronaut")
        return result

    def _batch_rows(self, iterator: Iterator[Dict[str, Any]]) -> Iterator[List[Dict[str, Any]]]:
        group_key = None
        group = []
        for row in iterator:
            key = [row[k_i] for k_i in self._key_columns]
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

    def create_index(self):
        key = None
        delimiter = self._delimiter.encode("ascii")
        with open_stream(self.filename, 'rb') as stream:
            header = stream.readline()
            header_cols = header.split(delimiter)
            column_keys = (
                header_cols.index(b'ModifiedPeptide'),
                header_cols.index(b'PrecursorCharge'),
            )
            offset = stream.tell()

            line = stream.readline()
            tokens = line.split(delimiter)
            key = (tokens[column_keys[0]].strip(b"_"), tokens[column_keys[1]])
            n = 0
            self.index.add(
                number=n,
                offset=offset,
                name=_ID_SEP.join(key).decode("utf8"),
                analyte=None
            )
            n += 1

            # To hold previous values
            last_offset = None
            while line:
                tokens = line.split(delimiter)
                next_key = (tokens[column_keys[0]].strip(b"_"), tokens[column_keys[1]])
                if next_key != key:
                    key = next_key
                    self.index.add(
                        number=n,
                        offset=offset,
                        name=_ID_SEP.join(key).decode("utf8"),
                        analyte=None
                    )
                    n += 1
                    last_offset = stream.tell()
                else:
                    offset = stream.tell()
                line = stream.readline()

        if key is not None:
            self.index.add(
                number=n,
                offset=last_offset,
                name=_ID_SEP.join(key).decode('utf8'),
                analyte=None
            )
        n += 1
        self.index.commit()
        return n

    def _generate_peaks(self, batch: List[Dict[str, Any]]) -> List[Tuple[float, float, List[annotation.IonAnnotationBase], List]]:
        peaks = []
        for row in batch:
            mz = float(row['FragmentMz'])
            intensity = float(row['RelativeIntensity'])

            series = row['FragmentType']
            ordinal = int(row['FragmentNumber'])
            charge = int(row['FragmentCharge'])

            loss_type = row['FragmentLossType']
            if loss_type != NO_LOSS:
                loss_type = ['-' + loss_type]
            else:
                loss_type = None

            annot = annotation.PeptideFragmentIonAnnotation(
                series, ordinal, neutral_losses=loss_type, charge=charge,
                mass_error=annotation.MassError(0, 'Da')
            )

            peak = [
                mz, intensity, [annot], []
            ]
            peaks.append(peak)
        return peaks

    def _build_analyte(self, description: Dict[str, Any], analyte: Analyte) -> Analyte:
        pf_seq = _rewrite_modified_peptide_as_proforma(description['ModifiedPeptide'])
        peptide = proforma.ProForma.parse(pf_seq)

        analyte.add_attribute(STRIPPED_PEPTIDE_TERM, description['StrippedPeptide'])
        analyte.add_attribute(PROFORMA_PEPTIDE_TERM, pf_seq)
        analyte.add_attribute("MS:1001117|theoretical mass", peptide.mass)
        analyte.add_attribute(CHARGE_STATE, int(description['PrecursorCharge']))

        protein_group_id = analyte.get_next_group_identifier()


        if 'UniProtIds' in description:
            analyte.add_attribute(
                "MS:1000885|protein accession",
                description['UniProtIds'],
                group_identifier=protein_group_id
            )
        if 'Protein Name' in description:
            analyte.add_attribute(
                "MS:1000886|protein name",
                description["Protein Name"],
                group_identifier=protein_group_id
            )
        if 'ProteinDescription' in description:
            analyte.add_attribute(
                "MS:1001088|protein description",
                description['ProteinDescription'],
                group_identifier=protein_group_id
            )

        if "OrganismId" in description:
            analyte.add_attribute_group([
                ["MS:1001467|taxonomy: NCBI TaxID", f"NCBITaxon:{description['OrganismId']}|{description['Organisms']}"],
                ["MS:1001469|taxonomy: scientific name", description['Organisms']],
            ])

        for key in self._custom_analyte_keys:
            if key in description:
                analyte.add_attribute_group([
                    [CUSTOM_ATTRIBUTE_NAME, key],
                    [CUSTOM_ATTRIBUTE_VALUE, _parse_value(description[key])]
                ])

        return analyte

    def _parse_from_buffer(self, buffer: List[Dict[str, Any]], spectrum_index: int = None) -> Spectrum:
        spec = self._new_spectrum()
        descr = buffer[0]

        key = (descr['ModifiedPeptide'].strip("_"), descr['PrecursorCharge'])

        spec.add_attribute(SPECTRUM_NAME, ID_SEP.join(key))
        spec.add_attribute(SELECTED_ION_MZ, float(descr['PrecursorMz']))
        # Charge does not belong in the spectrum
        # spec.add_attribute(CHARGE_STATE, int(descr['PrecursorCharge']))
        spec.add_attribute(SOURCE_FILE, descr['ReferenceRun'])
        spec.add_attribute(*self._spectrum_type())

        spec.add_attribute_group([
            [CUSTOM_ATTRIBUTE_NAME, "LabeledPeptide"],
            [CUSTOM_ATTRIBUTE_VALUE, descr['LabeledPeptide']]
        ])

        if 'IonMobility' in descr:
            spec.add_attribute("MS:1002476|ion mobility drift time", float(descr['IonMobility']))
        if 'CV' in descr:
            spec.add_attribute("MS:1001581|FAIMS compensation voltage", float(descr['CV']))
        if 'iRT' in descr:
            spec.add_attribute("MS:1000896|normalized retention time", float(descr['iRT']))

        analyte = self._new_analyte('1')
        self._build_analyte(descr, analyte)
        spec.add_analyte(analyte)

        for key in self._custom_spectrum_keys:
            if key in descr:
                spec.add_attribute_group([
                    [CUSTOM_ATTRIBUTE_NAME, key],
                    [CUSTOM_ATTRIBUTE_VALUE, _parse_value(descr[key])]
                ])

        spec.peak_list = self._generate_peaks(buffer)
        spec.add_attribute("MS:1003059|number of peaks", len(spec.peak_list))

        if spectrum_index is not None:
            spec.index = spectrum_index
        else:
            spec.index = -1
        spec.key = spec.index + 1
        return spec
