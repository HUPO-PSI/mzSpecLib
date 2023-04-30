import json

from typing import List, Tuple, Dict, Iterator, Any, Union

from pyteomics import proforma

from mzlib import annotation
from mzlib.backends.base import DEFAULT_VERSION, FORMAT_VERSION_TERM, _CSVSpectralLibraryBackendBase
from mzlib.spectrum import Spectrum, SPECTRUM_NAME


def _rewrite_unimod_peptide_as_proforma(sequence: str) -> str:
    return sequence.replace("(", '[').replace(')', ']').replace("UniMod", "UNIMOD")


CHARGE_STATE = "MS:1000041|charge state"
SELECTED_ION_MZ = "MS:1003208|experimental precursor monoisotopic m/z"
SOURCE_FILE = "MS:1003203|constituent spectrum file"
STRIPPED_PEPTIDE_TERM = "MS:1000888|stripped peptide sequence"
PROFORMA_PEPTIDE_TERM = "MS:1003169|proforma peptidoform sequence"

CUSTOM_ATTRIBUTE_NAME = "MS:1003275|other attribute name"
CUSTOM_ATTRIBUTE_VALUE = "MS:1003276|other attribute value"

NO_LOSS = 'noloss'


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


class DIANNTSVSpectralLibrary(_CSVSpectralLibraryBackendBase):
    format_name = "dia-nn.tsv"

    _custom_spectrum_keys = [
        "ExcludeFromAssay",
        "BGSInferenceId",
        "AllowForNormalization",
        "Workflow",
    ]

    _custom_analyte_keys = [
        "Proteotypic",
        "ProteinGroup",
    ]

    def __init__(self, filename: str, index_type=None, **kwargs):
        super().__init__(filename, index_type=index_type, delimiter='\t')

    def _spectrum_type(self):
        key = "MS:1003072|spectrum origin type"
        value = "MS:1003074|predicted spectrum"
        return key, value

    def read_header(self) -> bool:
        result = super().read_header()
        self.add_attribute(FORMAT_VERSION_TERM, DEFAULT_VERSION)
        self.add_attribute("MS:1003207|library creation software", "MS:1003253|DIA-NN")
        return result

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

    def _parse_from_buffer(self, buffer: List[Dict[str, Any]], spectrum_index: int = None) -> Spectrum:
        spec = self._new_spectrum()
        descr = buffer[0]

        spec.add_attribute(SPECTRUM_NAME, descr['transition_group_id'])
        spec.add_attribute(SELECTED_ION_MZ, float(descr['PrecursorMz']))
        spec.add_attribute(CHARGE_STATE, int(descr['PrecursorCharge']))
        spec.add_attribute(SOURCE_FILE, descr['FileName'])
        spec.add_attribute(*self._spectrum_type())

        if int(descr['decoy']):
            spec.add_attribute("MS:1003072|spectrum origin type", "MS:1003192|decoy spectrum")

        if 'IonMobility' in descr:
            spec.add_attribute("MS:1002476|ion mobility drift time", float(descr['IonMobility']))

        analyte = self._new_analyte('1')

        pf_seq = _rewrite_unimod_peptide_as_proforma(descr['FullUniModPeptideName'])
        peptide = proforma.ProForma.parse(pf_seq)

        analyte.add_attribute(STRIPPED_PEPTIDE_TERM, descr['PeptideSequence'])
        analyte.add_attribute(PROFORMA_PEPTIDE_TERM, pf_seq)
        analyte.add_attribute("MS:1001117|theoretical mass", peptide.mass)

        protein_group_id = analyte.get_next_group_identifier()
        if "UniprotID" in descr:
            analyte.add_attribute(
                "MS:1000885|protein accession",
                descr['UniprotID'],
                group_identifier=protein_group_id
            )
        if "ProteinName" in  descr:
            analyte.add_attribute(
                "MS:1000886|protein name",
                descr["Protein Name"],
                group_identifier=protein_group_id
            )

        for key in self._custom_analyte_keys:
            if key in descr:
                analyte.add_attribute_group([
                    [CUSTOM_ATTRIBUTE_NAME, key],
                    [CUSTOM_ATTRIBUTE_VALUE, _parse_value(descr[key])]
                ])

        spec.add_analyte(analyte)
        spec.peak_list = self._generate_peaks(buffer)
        spec.add_attribute("MS:1003059|number of peaks", len(spec.peak_list))

        for key in self._custom_spectrum_keys:
            if key in descr:
                spec.add_attribute_group([
                    [CUSTOM_ATTRIBUTE_NAME, key],
                    [CUSTOM_ATTRIBUTE_VALUE, _parse_value(descr[key])]
                ])

        if spectrum_index is not None:
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

    def _batch_rows(self, iterator: Iterator[Dict[str, Any]]) -> Iterator[List[Dict[str, Any]]]:
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


DiaNNTSVSpectralLibrary = DIANNTSVSpectralLibrary