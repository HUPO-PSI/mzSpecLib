import re
import io
import os
import logging

from typing import Any, Callable, Collection, Dict, List, Mapping, Optional, Set, Tuple, Iterable
import warnings

from pyteomics import proforma

from mzlib import annotation

from mzlib.analyte import FIRST_ANALYTE_KEY, FIRST_INTERPRETATION_KEY, Analyte
from mzlib.spectrum import Spectrum, SPECTRUM_NAME
from mzlib.attributes import AttributeManager, Attributed

from .base import DEFAULT_VERSION, FORMAT_VERSION_TERM, _PlainTextSpectralLibraryBackendBase, LIBRARY_NAME_TERM
from .utils import try_cast, open_stream


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


leader_terms = {
    "Name": SPECTRUM_NAME,
}

STRIPPED_PEPTIDE_TERM = "MS:1000888|stripped peptide sequence"

analyte_terms = {
    "MW": "MS:1000224|molecular mass",
    "ExactMass": "MS:1000224|molecular mass",
    "Theo_mz_diff": "MS:1003209|monoisotopic m/z deviation",
    "Scan": {
        "Protein": "MS:1000885|protein accession",
        "Mods": "MS:1001471|peptide modification details",
        "Naa": "MS:1003043|number of residues",
    },
    "Pep": {
        "Tryptic": [["MS:1003048|number of enzymatic termini", 2], ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
        "SemiTryptic": [["MS:1003048|number of enzymatic termini", 1], ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
        "N-Semitryptic": [["MS:1003048|number of enzymatic termini", 1], ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
        "C-Semitryptic": [["MS:1003048|number of enzymatic termini", 1], ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
        "Tryptic/miss_good_confirmed": [["MS:1003048|number of enzymatic termini", 2],
                                        ["MS:1003044|number of missed cleavages", "0"],
                                        ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
        "Tryptic/miss_bad_confirmed": [["MS:1003048|number of enzymatic termini", 2],
                                        ["MS:1003044|number of missed cleavages", ">0"],
                                        ["MS:1001045|cleavage agent name", "MS:1001251|Trypsin"]],
        },
    "MC": "MS:1003044|number of missed cleavages",
    # "Protein": "MS:1000885|protein accession",
    "Mods": "MS:1001471|peptide modification details",
    "Naa": "MS:1003043|number of residues",
    "PrecursorMonoisoMZ": "MS:1003208|experimental precursor monoisotopic m/z",
    "Mz_exact": "MS:1003208|experimental precursor monoisotopic m/z",
    "Mz_av": "MS:1003054|theoretical average m/z",
}


other_terms = {
    "Charge": "MS:1000041|charge state",
    "Parent": "MS:1000744|selected ion m/z",
    "ObservedPrecursorMZ": "MS:1000744|selected ion m/z",
    "Single": ["MS:1003065|spectrum aggregation type", "MS:1003066|singleton spectrum"],
    "Consensus": ["MS:1003065|spectrum aggregation type", "MS:1003067|consensus spectrum"],
    "Inst": {"it": [["MS:1000044|dissociation method", "MS:1002472|trap-type collision-induced dissociation"]],
             "hcd": [["MS:1000044|dissociation method", "MS:1000422|beam-type collision-induced dissociation"]]},
    "Spec": {"Consensus": [["MS:1003065|spectrum aggregation type", "MS:1003067|consensus spectrum"]]},
    "Scan": "MS:1003057|scan number",
            "Origfile": "MS:1003203|constituent spectrum file",
            "Sample": "MS:1000002|sample name",
            "Filter": "MS:1000512|filter string",
            "FTResolution": "MS:1000028|detector resolution",
            "ms1PrecursorAb": "MS:1003085|previous MS1 scan precursor intensity",
            "Precursor1MaxAb": "MS:1003086|precursor apex intensity",
            "Purity": "MS:1009013|isolation window precursor purity",
            "Unassigned": "MS:1003080|top 20 peak unassigned intensity fraction",
            "Unassign_all": "MS:1003079|total unassigned intensity fraction",
            "BasePeak": "MS:1000505|base peak intensity",
            "Num peaks": "MS:1003059|number of peaks",
}


interpretation_member_terms = {
    "Q-value": "MS:1002354|PSM-level q-value"
}


species_map = {
    "human": [["MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:9606|Homo sapiens"],
                ["MS:1001469|taxonomy: scientific name",
                    "Homo sapiens"],
                ["MS:1001468|taxonomy: common name", "human"]],
    "zebrafish": [["MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:7955|Danio rerio"],
                    ["MS:1001469|taxonomy: scientific name",
                        "Danio rerio"],
                    ["MS:1001468|taxonomy: common name", "zebra fish"]],
    "chicken": [["MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:9031|Gallus gallus"],
                ["MS:1001469|taxonomy: scientific name",
                    "Gallus gallus"],
                ["MS:1001468|taxonomy: common name", "chicken"]],
}


MODIFICATION_NAME_MAP = {
    "CAM": "Carbamidomethyl",
    "Pyro_glu": "Glu->pyro-Glu", # Resolves UNIMOD ambiguity
    "Pyro-glu": "Gln->pyro-Glu",
    "Oxidation": "Oxidation",
    "Phospho": "Phospho",
    "TMT6plex": "TMT6plex",
    "iTRAQ": "iTRAQ",
    "Acetyl": "Acetyl",
    "TMT": "TMT6plex",
}


# TODO: ppm is unsigned, add mass calculation to determine true mass accuracy

annotation_pattern = re.compile(r"""^
(?:(?:(?P<series>[abyxcz]\.?)(?P<ordinal>\d+))|
   (:?Int/(?P<series_internal>[ARNDCEQGHKMFPSTWYVILJarndceqghkmfpstwyvilj]+))|
   (?P<precursor>p)|
   (:?I(?P<immonium>[ARNDCEQGHKMFPSTWYVIL])(?:(?P<immonium_modification>CAM)|[A-Z])?)|
   (?P<reporter>(?:TMT|iTRAQ)(?P<reporter_mass>\d+[NC]?))|
   (?:_(?P<external_ion>[^\s,/]+))|
   (?P<unannotated>\?)
)
(?P<neutral_losses>(?:[+-]\d*[A-Z][A-Za-z0-9]*)+)?
(?:\[M(?P<adducts>(:?[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?
(?:\^(?P<charge>[+-]?\d+))?
(?:(?P<isotope>[+-]\d*)i)?
(?:@(?P<analyte_reference>[^/\s]+))?
(?:/(?P<mass_error>[+-]?\d+(?:\.\d+))(?P<mass_error_unit>ppm)?)?
""", re.X)


class MSPAnnotationStringParser(annotation.AnnotationStringParser):
    def _dispatch_internal_peptide_fragment(self, data: Dict[str, Any], adducts: List, charge: int, isotope: int, neutral_losses: List,
                                            analyte_reference: Any, mass_error: Any, **kwargs):
        spectrum = kwargs.get("spectrum")
        if spectrum is None:
            raise ValueError("Cannot infer sequence coordinates from MSP internal fragmentation notation without"
                             " a reference to the spectrum, please pass spectrum as a keyword argument")
        sequence = self._get_peptide_sequence_for_analyte(
            spectrum, analyte_reference)
        subseq = data['series_internal'].upper()
        try:
            start_index = sequence.index(subseq)
        except ValueError as err:
            raise ValueError(
                f"Cannot locate internal subsequence {subseq} in {sequence}") from err
        end_index = start_index + len(subseq)
        data['internal_start'] = start_index + 1
        data['internal_end'] = end_index
        return super(MSPAnnotationStringParser, self)._dispatch_internal_peptide_fragment(
            data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, **kwargs)

    def _dispatch_immonium(self, data: Dict[str, Any], adducts: List, charge: int, isotope: int, neutral_losses: List,
                           analyte_reference: Any, mass_error: Any, **kwargs):
        modification = data['immonium_modification']
        if modification is not None:
            try:
                modification = MODIFICATION_NAME_MAP[modification]
                data['immonium_modification'] = modification
            except KeyError as err:
                print(f"Failed to convert immonium ion modification {modification}")
        return super(MSPAnnotationStringParser, self)._dispatch_immonium(
            data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, **kwargs)

    def _get_peptide_sequence_for_analyte(self, spectrum: Spectrum, analyte_reference: Optional[Any]=None) -> str:
        if analyte_reference is None:
            if len(spectrum.analytes) == 0:
                return None
            else:
                analyte_reference = spectrum.analytes[FIRST_ANALYTE_KEY].id
        analyte = spectrum.analytes.get(analyte_reference)
        if analyte is None:
            return None
        return analyte.get_attribute('MS:1000888|stripped peptide sequence')


parse_annotation = MSPAnnotationStringParser(annotation_pattern)
MODIFICATION_LIST_PARSER = re.compile(r"(\d+),([ARNDCEQGHKMFPSTWYVIL_\-]),([A-Za-z0-9_\-]+)")


class ModificationParser:
    pattern: re.Pattern
    modification_map: Dict[str, str]
    unknown_modifications: Set[str]

    def __init__(self, pattern: str=None, modification_map: Dict[str, str]=None):
        if pattern is None:
            pattern = MODIFICATION_LIST_PARSER.pattern
        if modification_map is None:
            modification_map = dict(MODIFICATION_NAME_MAP)
        self.pattern = re.compile(pattern)
        self.modification_map = modification_map
        self.unknown_modifications = set()

    def __call__(self, text: str):
        return self.parse(text)

    def parse(self, text: str) -> List[Tuple[int, str, str]]:
        if not isinstance(text, str) or not text:
            return []
        i = 0
        n = len(text)

        mods = []
        while i < n and text[i].isdigit():
            i += 1

        for position, residue, mod in self.pattern.findall(text):
            position = int(position)
            if mod not in self.modification_map:
                warnings.warn(f"{mod} is not found in the known MSP modification mapping. Using this name verbatim")
                modification_name = mod
                self.modification_map[mod] = mod
                self.unknown_modifications.add(mod)
            else:
                modification_name = self.modification_map[mod]
            mods.append((position, residue, modification_name))
        return mods


class AttributeHandler:
    keys: Collection[str]

    def __init__(self, keys: Collection[str]):
        self.keys = keys

    def __contains__(self, key):
        return key in self.keys

    def add_value(self, key: str, value: Any, container: Attributed):
        container.add_attribute(key, value)

    def add_group(self, keys: List[str], values: List[Any], container: Attributed):
        group_id = container.get_next_group_identifier()
        for k, v in zip(keys, values):
            container.add_attribute(k, v, group_id)

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        raise NotImplementedError()

    def __call__(self, key: str, value: Any, container: Attributed) -> bool:
        return self.handle(key, value, container)

    def chain(self, handler: 'AttributeHandler') -> 'AttributeHandlerChain':
        return AttributeHandlerChain([self, handler])

    def __and__(self, handler: 'AttributeHandler') -> 'AttributeHandlerChain':
        return self.chain(handler)


class MappingAttributeHandler(AttributeHandler):
    keys: Dict[str, Any]

    def __init__(self, keys: Dict[str, Any]):
        self.keys = keys

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        trans_key = self.keys[key]
        if value is None:
            if isinstance(trans_key, list):
                k, v = trans_key
                v = try_cast(v)
                self.add_value(k, v, container)
            else:
                return False
        elif isinstance(trans_key, str):
            self.add_value(trans_key, value, container)
        elif isinstance(trans_key, dict):
            if value in trans_key:
                # If the mapping is a plain string, add it
                if isinstance(trans_key[value], str):
                    key, value = trans_key[value].split("=")
                    self.add_value(key, try_cast(value), container)
                # Or if it is a list, then there are multiple terms to add within a group
                elif isinstance(trans_key[value], list):
                    if len(trans_key[value]) == 1:
                        for item in trans_key[value]:
                            self.add_value(item[0], try_cast(item[1]), container)
                    else:
                        k, v = zip(*trans_key[value])
                        self.add_group(k, v, container)
            else:
                return False
        return True


class AttributeHandlerChain:
    chain: List[AttributeHandler]

    def __init__(self, chain: List[AttributeHandler]):
        self.chain = chain

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        result = True
        for handler in self.chain:
            result |= handler.handle(key, value, container)
            if not result:
                return result
        return result

    def __contains__(self, key):
        for handler in self.chain:
            if key in handler:
                return True
        return False

    def __iter__(self):
        return iter(self.chain)

    def __call__(self, key: str, value: Any, container: Attributed) -> bool:
        return self.handle(key, value, container)

    def chain(self, handler: 'AttributeHandler') -> 'AttributeHandlerChain':
        return self.__class__(self.chain + [handler])

    def add(self, handler: AttributeHandler):
        self.chain.append(handler)
        return handler


class RegexAttributeHandler(AttributeHandler):
    pattern: re.Pattern
    variable_attributes: List[str]
    constant_attributes: Optional[List[Tuple[str, Any]]]

    def __init__(self, keys: Collection[str], pattern: re.Pattern, variable_attributes, constant_attributes=None):
        super().__init__(keys)
        self.pattern = pattern
        self.variable_attributes = variable_attributes
        self.constant_attributes = constant_attributes

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        match = self.pattern.match(value)
        if match:
            if self.constant_attributes:
                ks, vs = map(list, zip(*self.constant_attributes))
                ks.extend(self.variable_attributes)
                vs.extend(match.groups())
                self.add_group(ks, vs, container)
            else:
                self.add_value(self.variable_attributes[0], match.group(1), container)
            return True
        return False


class FunctionAttributeHandler(AttributeHandler):
    func: Callable[[str, Any, Attributed], bool]

    def __init__(self, keys: Collection[str], func: Callable[[str, Any, Attributed], bool]):
        super().__init__(keys)
        self.func = func

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        return self.func(key, value, container)

    @classmethod
    def wraps(cls, *keys):
        def wrapper(func):
            return cls(keys, func)
        return wrapper


class DispatchingAttributeHandler(AttributeHandlerChain):
    mapping: Dict[str, AttributeHandler]

    def __init__(self, chain: List[AttributeHandler]=None):
        if not chain:
            chain = []
        super().__init__(chain)
        self.mapping = {}
        for handler in self:
            for key in handler.keys:
                self.mapping[key] = handler

    def handle(self, key: str, value: Any, container: Attributed) -> bool:
        handler = self.mapping[key]
        return handler(key, value, container)

    def __contains__(self, key):
        return key in self.mapping

    def add(self, handler: AttributeHandler):
        super().add(handler)
        for key in handler.keys:
            self.mapping[key] = handler
        return handler


msp_spectrum_attribute_handler = DispatchingAttributeHandler()
msp_analyte_attribute_handler = DispatchingAttributeHandler()

@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("HCD")
def dissociation_method_handler(key: str, value: str, container: Attributed) -> bool:
    container.add_attribute("MS:1000044|dissociation method",
                            "MS:1000422|beam-type collision-induced dissociation")
    found_match = False
    if value is not None:
        match = re.match(r"([\d\.]+)\s*ev", value, flags=re.IGNORECASE)
        if match is not None:
            found_match = True
            group_identifier = container.get_next_group_identifier()
            container.add_attribute("MS:1000045|collision energy",
                                try_cast(match.group(1)), group_identifier)
            container.add_attribute("UO:0000000|unit", "UO:0000266|electronvolt", group_identifier)
        match = re.match(r"([\d\.]+)\s*%", value)
        if match is not None:
            found_match = True
            group_identifier = container.get_next_group_identifier()
            container.add_attribute("MS:1000045|collision energy",
                                try_cast(match.group(1)), group_identifier)
            container.add_attribute("UO:0000000|unit", "UO:0000187|percent", group_identifier)
    return found_match


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("Collision_energy", "CE")
def collision_energy_handler(key: str, value: str, container: Attributed) -> bool:
    if isinstance(value, str):
        match = re.match(r"([\d\.]+)", value)
        if match is not None:
            value = try_cast(match.group(1))
    if value is not None:
        if match is not None:
            group_identifier = container.get_next_group_identifier()
            container.add_attribute(
                "MS:1000045|collision energy", value, group_identifier)
            container.add_attribute(
                "UO:0000000|unit", "UO:0000266|electronvolt", group_identifier)
            return True
    return False


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("RT")
def rt_handler(key, value, container) -> bool:
    match = re.match(r"([\d\.]+)\s*(\D*)",
                     value)
    if match is not None:
        if match.group(2):
            container.add_attribute(
                "ERROR", f"Need more RT parsing code to handle this value")
            return False
        else:
            group_identifier = container.get_next_group_identifier()
            container.add_attribute(
                "MS:1000894|retention time", try_cast(match.group(1)), group_identifier)
            #### If the value is greater than 250, assume it must be seconds
            if float(match.group(1)) > 250:
                container.add_attribute(
                    "UO:0000000|unit", "UO:0000010|second", group_identifier)
            #### Although normally assume minutes
            else:
                container.add_attribute(
                    "UO:0000000|unit", "UO:0000031|minute", group_identifier)
            return True
    return False


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("ms2IsolationWidth")
def isolation_width_handler(key, value, container) -> bool:
    if value is None:
        return False
    group_identifier = container.get_next_group_identifier()
    container.add_attribute("MS:1000828|isolation window lower offset",
        (float(value) / 2), group_identifier)
    container.add_attribute("UO:0000000|unit",
                        "MS:1000040|m/z", group_identifier)
    group_identifier = container.get_next_group_identifier()
    container.add_attribute("MS:1000829|isolation window upper offset",
        (float(value) / 2), group_identifier)
    container.add_attribute("UO:0000000|unit",
                        "MS:1000040|m/z", group_identifier)
    return True


@msp_analyte_attribute_handler.add
@FunctionAttributeHandler.wraps("Mz_diff", "Theo_mz_diff")
def mz_diff_handler(key, value, container: Attributed) -> bool:
    if isinstance(value, float):
        # We must be dealing with a unit-less entry.
        group_identifier = container.get_next_group_identifier()
        container.add_attribute(
            "MS:1001975|delta m/z", value, group_identifier)
        container.add_attribute(
            "UO:0000000|unit", "MS:1000040|m/z", group_identifier)
    else:
        match = re.match(
            r"([\-\+e\d\.]+)\s*ppm", value, flags=re.IGNORECASE)
        if match is not None:
            group_identifier = container.get_next_group_identifier()
            container.add_attribute(
                "MS:1001975|delta m/z", try_cast(match.group(1)), group_identifier)
            container.add_attribute(
                "UO:0000000|unit", "UO:0000169|parts per million", group_identifier)
        else:
            match = re.match(
                r"([\-\+e\d\.]+)\s*", value)
            if match is not None:
                group_identifier = container.get_next_group_identifier()
                container.add_attribute(
                    "MS:1001975|delta m/z", try_cast(match.group(1)), group_identifier)
                container.add_attribute(
                    "UO:0000000|unit", "MS:1000040|m/z", group_identifier)
            else:
                return False
    return True


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("Dev_ppm")
def dev_ppm_handler(key, value, container) -> bool:
    if value is None:
        return False
    group_identifier = container.get_next_group_identifier()
    container.add_attribute(
        "MS:1001975|delta m/z", try_cast(value), group_identifier)
    container.add_attribute(
        "UO:0000000|unit", "UO:0000169|parts per million", group_identifier)
    return True


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("Nrep", "Nreps")
def nreps_handler(key, value, container):
    if value is None:
        return False
    if not isinstance(value, str):
        value = str(value)
    match = re.match(
        r"(\d+)/(\d+)", value)
    if match is not None:
        container.add_attribute(
            "MS:1003070|number of replicate spectra used", try_cast(match.group(1)))
        container.add_attribute(
            "MS:1003069|number of replicate spectra available", try_cast(match.group(2)))
        return True
    else:
        match = re.match(
            r"(\d+)", value)
        if match is not None:
            container.add_attribute(
                "MS:1003070|number of replicate spectra used", try_cast(match.group(1)))
            container.add_attribute(
                "MS:1003069|number of replicate spectra available", try_cast(match.group(1)))
            return True
        else:
            return False


@msp_analyte_attribute_handler.add
@FunctionAttributeHandler.wraps("Organism")
def organism_handler(key, value, container):
    if value is None:
        return False
    value = value.strip('"')

    if value in species_map:
        group_identifier = container.get_next_group_identifier()
        for item in species_map[value]:
            container.add_attribute(
                item[0], try_cast(item[1]), group_identifier)
        return True
    return False


@msp_analyte_attribute_handler.add
@FunctionAttributeHandler.wraps("Protein")
def protein_handler(key, value, container):
    if value is None:
        return False
    key = "MS:1000885|protein accession"
    match = re.match(r"\(pre=(.),post=(.)\)", value)
    if match is not None:
        value = value[:match.start()]
        container.add_attribute("MS:1001112|n-terminal flanking residue", match.group(1))
        container.add_attribute("MS:1001113|c-terminal flanking residue", match.group(2))
    container.add_attribute(key, value)
    return True


class MSPSpectralLibrary(_PlainTextSpectralLibraryBackendBase):
    file_format = "msp"
    format_name = "msp"

    modification_parser: ModificationParser
    unknown_attributes: Set[str]

    def __init__(self, filename, index_type=None, read_metadata=True):
        super().__init__(filename, index_type, read_metadata)
        self.modification_parser = ModificationParser()
        self.unknown_attributes = set()

    @classmethod
    def guess_from_header(cls, filename: str) -> bool:
        with open_stream(filename, 'r') as stream:
            first_line = stream.readline()
            if re.match("Name: ", first_line):
                return True
        return False

    def read_header(self) -> bool:
        with open_stream(self.filename, 'r') as stream:
            match, offset = self._parse_header_from_stream(stream)
            return match
        return False

    def _parse_header_from_stream(self, stream: io.IOBase) -> Tuple[bool, int]:
        first_line = stream.readline()
        attributes = AttributeManager()
        attributes.add_attribute(FORMAT_VERSION_TERM, DEFAULT_VERSION)
        attributes.add_attribute(LIBRARY_NAME_TERM, self.filename)
        self.attributes.clear()
        self.attributes._from_iterable(attributes)
        if re.match("Name: ", first_line):
            return True, 0
        return False, 0

    def create_index(self) -> int:
        """
        Populate the spectrum index

        Returns
        -------
        n_spectra: int
            The number of entries read
        """

        #### Check that the spectrum library filename isvalid
        filename = self.filename

        #### Determine the filesize
        file_size = os.path.getsize(filename)

        with open_stream(filename, 'r') as infile:
            state = 'header'
            spectrum_buffer = []
            n_spectra = 0
            start_index = 0
            file_offset = 0
            line_beginning_file_offset = 0
            spectrum_file_offset = 0
            spectrum_name = ''

            # Required for counting file_offset manually (LF vs CRLF)
            infile.readline()
            file_offset_line_ending = len(infile.newlines) - 1
            infile.seek(0)

            # if debug:
            #     eprint("INFO: Reading..", end='', flush=True)
            logger.debug(f"Reading {filename} ({file_size} bytes)...")
            while 1:
                line = infile.readline()
                if len(line) == 0:
                    break

                line_beginning_file_offset = file_offset

                #### tell() is twice as slow as counting it myself
                # file_offset = infile.tell()
                file_offset += len(line) + file_offset_line_ending

                line = line.rstrip()
                if state == 'header':
                    if re.match('Name: ', line):
                        state = 'body'
                        spectrum_file_offset = line_beginning_file_offset
                    else:
                        continue
                if state == 'body':
                    if len(line) == 0:
                        continue
                    if re.match('Name: ', line):
                        if len(spectrum_buffer) > 0:
                            self.index.add(
                                number=n_spectra + start_index,
                                offset=spectrum_file_offset,
                                name=spectrum_name,
                                analyte=None)
                            n_spectra += 1
                            spectrum_buffer = []
                            #### Commit every now and then
                            if n_spectra % 10000 == 0:
                                self.index.commit()
                                logger.info(f"... Indexed  {file_offset} bytes, {n_spectra} spectra read")

                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = re.match(r'Name:\s+(.+)', line).group(1)

                    spectrum_buffer.append(line)
            logger.debug(f"Processed {file_offset} bytes, {n_spectra} spectra read")
            self.index.add(
                number=n_spectra + start_index,
                offset=spectrum_file_offset,
                name=spectrum_name,
                analyte=None)
            self.index.commit()
            n_spectra += 1

            #### Flush the index
            self.index.commit()

        return n_spectra

    def _buffer_from_stream(self, infile: Iterable[str]) -> List[str]:
        state = 'body'
        spectrum_buffer = []

        file_offset = 0
        line_beginning_file_offset = 0
        spectrum_file_offset = 0
        spectrum_name = ''
        for line in infile:
            line_beginning_file_offset = file_offset
            file_offset += len(line)
            line = line.rstrip()
            if state == 'body':
                if len(line) == 0:
                    continue
                if re.match(r'Name: ', line):
                    if len(spectrum_buffer) > 0:
                        return spectrum_buffer
                    spectrum_file_offset = line_beginning_file_offset
                    spectrum_name = re.match(r'Name:\s+(.+)', line).group(1)
                spectrum_buffer.append(line)
        return spectrum_buffer

    def _parse(self, buffer: Iterable[str], spectrum_index: int=None) -> Spectrum:

        #### Start in the header section of the entry
        in_header = True

        #### Reset all spectrum properties in case there is already data here
        attributes = {}
        peak_list = []

        #### Loop through each line in the buffered list
        for line in buffer:

            #print(line)

            #### If in the the header portion of the entry
            if in_header:

                #### Extract the key,value pair by splitting on the *first* colon with optional whitespace
                match = re.match(r"\s*#", line)
                if match:
                    continue
                elif line.count(":") > 0:
                    key, value = re.split(r":\s*", line, 1)
                elif line.count("=") > 0:
                    key, value = re.split(r"=\s*", line, 1)
                elif line.count("\t") > 0:
                    logger.error("Looks like peaks in the header???")
                    in_header = False
                else:
                    key = line
                    value = None

                #### Adds the key-value pair to the dictionary
                attributes[key] = value

                #### If the key is "Num peaks" then we're done with the header and peak list follows
                if key == "Num peaks":
                    in_header = False

                #### The "Comment" key requires special parsing
                if key == "Comment":

                    #### Remove it from attributes
                    del attributes[key]
                    self._parse_comment(value, attributes)

            #### Else in the peaks section. Parse the peaks.
            else:
                #### Split into the expected three values
                values = re.split(r'\s+', line, maxsplit=2)
                interpretations = ""
                aggregation = ""
                if len(values) == 1:
                    mz = values
                    intensity = "1"
                if len(values) == 2:
                    mz, intensity = values
                elif len(values) == 3:
                    mz, intensity, interpretations = values
                elif len(values) > 3:
                    mz, intensity, interpretations = values[0:2]
                else:
                    mz = "1"
                    intensity = "1"

                interpretations = interpretations.strip('"')
                if interpretations.startswith("?"):
                    parts = re.split(r"\s+", interpretations)
                    if len(parts) > 1:
                        # Some msp files have a concept for ?i, but this requires a definition
                        interpretations = "?"
                        aggregation = parts[1:]
                else:
                    if " " in interpretations:
                        parts = re.split(r"\s+", interpretations)
                        interpretations = parts[0]
                        aggregation = parts[1:]

                #### Add to the peak list
                peak_list.append([float(mz), float(intensity), interpretations, aggregation])

        #### Now convert the format attributes to standard ones
        spectrum = self._make_spectrum(peak_list, attributes)
        if spectrum_index is not None:
            spectrum.index = spectrum_index
        else:
            spectrum.index = -1
        spectrum.key = spectrum.index + 1
        for i, peak in enumerate(spectrum.peak_list):
            try:
                parsed_interpretation = self._parse_annotation(peak[2], spectrum=spectrum)
            except ValueError as err:
                message = str(err)
                raise ValueError(
                    f"An error occurred while parsing the peak annotation for peak {i}: {message}") from err
            peak[2] = parsed_interpretation

        return spectrum

    def _parse_annotation(self, annotation: str, wrap_errors: bool=True, **kwargs):
        return parse_annotation(annotation_string=annotation, wrap_errors=wrap_errors, **kwargs)

    def _parse_comment(self, value: str, attributes: Attributed):
        comment_items = re.split(" ", value)

        #### Any spaces within quotes are then de-split
        fixed_comment_items = []
        quote_counter = 0
        new_item = ""
        for item in comment_items:
            if new_item > "":
                new_item = new_item + " "
            new_item = new_item + item
            n_quotes = new_item.count('"')
            if n_quotes/2 == int(n_quotes/2):
                fixed_comment_items.append(new_item)
                new_item = ""

        #### Try to split each item on the first = character and store in attributes
        for item in fixed_comment_items:
            #### If there is an = then split on the first one and store key and value
            if item.count("=") > 0:
                comment_key, comment_value = item.split("=", 1)
                attributes[comment_key] = try_cast(comment_value)

            #### Otherwise just store the key with a null value
            else:
                attributes[item] = None

    def _make_attribute_handlers(self):
        other_manager = MappingAttributeHandler(other_terms)
        analyte_manager = MappingAttributeHandler(analyte_terms)
        interpretation_member_manager = MappingAttributeHandler(interpretation_member_terms)
        return (other_manager,
                analyte_manager,
                interpretation_member_manager,
                msp_spectrum_attribute_handler,
                msp_analyte_attribute_handler)

    def _make_spectrum(self, peak_list: List, attributes: Mapping[str, str]):
        spectrum = self._new_spectrum()
        interpretation = self._new_interpretation(FIRST_INTERPRETATION_KEY)
        interpretation_member = None
        analyte = self._new_analyte(FIRST_ANALYTE_KEY)
        spectrum.add_analyte(analyte)
        # interpretation.add_analyte(analyte)
        spectrum.peak_list = peak_list
        # spectrum.interpretations.add_interpretation(interpretation)

        #### Add special terms that we want to start off with
        for term in leader_terms:
            if term in attributes:
                spectrum.add_attribute(
                    leader_terms[term], try_cast(attributes[term]))
            else:
                spectrum.add_attribute(
                    "ERROR", f"Required term {leader_terms[term]} is missing")

        #### Translate the rest of the known attributes and collect unknown ones
        unknown_terms = []

        (other_manager, analyte_manager, interpretation_member_manager,
         spectrum_func_handler, analyte_func_handler) = self._make_attribute_handlers()
        for attribute in attributes:

            #### Skip a leader term that we already processed
            if attribute in leader_terms:
                continue

            if attribute in other_manager:
                if not other_manager(attribute, attributes[attribute], spectrum):
                    unknown_terms.append(attribute)
            elif attribute in analyte_manager:
                if not analyte_manager(attribute, attributes[attribute], analyte):
                    unknown_terms.append(attribute)
            elif attribute in interpretation_member_manager:
                if interpretation_member is None:
                    interpretation_member = self._new_interpretation_member(1)
                    interpretation.add_member_interpretation(interpretation_member)
                if not interpretation_member_manager(attribute, attributes[attribute], interpretation_member):
                    unknown_terms.append(attribute)

            elif attribute in spectrum_func_handler:
                if not spectrum_func_handler(attribute, attributes[attribute], spectrum):
                    unknown_terms.append(attribute)

            elif attribute in analyte_func_handler:
                if not analyte_func_handler(attribute, attributes[attribute], analyte):
                    unknown_terms.append(attribute)

            #### Expand the Fullname attribute
            elif attribute == "Fullname":
                if attributes[attribute] is not None:
                    match = re.match(
                        r"([A-Z\-\*])\.([A-Z]+)\.([A-Z\-\*])/*([\d]*)", attributes[attribute])
                    if match is not None:
                        analyte.add_attribute(
                            "MS:1000888|stripped peptide sequence", match.group(2))
                        analyte.add_attribute(
                            "MS:1001112|n-terminal flanking residue", match.group(1))
                        analyte.add_attribute(
                            "MS:1001113|c-terminal flanking residue", match.group(3))
                        if match.group(4):
                            spectrum.add_attribute(
                                "MS:1000041|charge state", try_cast(match.group(4)))
                    else:
                        spectrum.add_attribute(
                            "ERROR", f"Unable to parse {attributes[attribute]} in {attribute}")
                        unknown_terms.append(attribute)
                else:
                    spectrum.add_attribute(
                        "ERROR", f"Attribute {attribute} must have a value")
                    unknown_terms.append(attribute)

            #### Otherwise add this term to the list of attributes that we don't know how to handle
            else:
                unknown_terms.append(attribute)

        if not analyte.has_attribute(STRIPPED_PEPTIDE_TERM):
            if spectrum.has_attribute(SPECTRUM_NAME):
                name = spectrum.get_attribute(SPECTRUM_NAME)
                # lookup = spectrum.attribute_dict[SPECTRUM_NAME]
                # name = spectrum.attributes[lookup["indexes"][0]][1]
                match = re.match(r"(.+)/(\d+)", name)
                if match:
                    analyte.add_attribute(
                        STRIPPED_PEPTIDE_TERM, match.group(1))
                    spectrum.add_attribute(
                        "MS:1000041|charge state", try_cast(match.group(2)))

        #### Handle the uninterpretable terms
        for attribute in unknown_terms:
            self.unknown_attributes.add(attribute)
            if attributes[attribute] is None:
                spectrum.add_attribute(
                    "MS:1003275|other attribute name", try_cast(attribute))
            else:
                group_identifier = spectrum.get_next_group_identifier()
                spectrum.add_attribute(
                    "MS:1003275|other attribute name", try_cast(attribute), group_identifier)
                spectrum.add_attribute("MS:1003276|other attribute value",
                                   try_cast(attributes[attribute]), group_identifier)

        self._complete_analyte(analyte)

        if analyte:
            spectrum.add_analyte(analyte)
            interpretation.add_analyte(analyte)
            spectrum.add_interpretation(interpretation)
        return spectrum

    def _complete_analyte(self, analyte: Analyte):
        peptide = proforma.ProForma.parse(analyte.get_attribute("MS:1000888|stripped peptide sequence"))
        if analyte.has_attribute("MS:1001471|peptide modification details"):
            modification_details = analyte.get_attribute("MS:1001471|peptide modification details")
            mods = self.modification_parser(modification_details)
            for position, residue, mod in mods:
                seqpos = list(peptide.sequence[position])
                if not seqpos[1]:
                    seqpos[1] = [proforma.GenericModification(mod)]
                peptide.sequence[position] = tuple(seqpos)
                assert seqpos[0] == residue
        analyte.add_attribute("MS:1003169|proforma peptidoform sequence", str(peptide))
        analyte.remove_attribute("MS:1001471|peptide modification details")
        analyte.add_attribute("MS:1001117|theoretical mass", peptide.mass)

    def get_spectrum(self, spectrum_number: int=None, spectrum_name: str=None) -> Spectrum:
        # keep the two branches separate for the possibility that this is not possible with all
        # index schemes.
        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError("Provide only one of spectrum_number or spectrum_name")
            offset = self.index.offset_for(spectrum_number)
        elif spectrum_name is not None:
            offset = self.index.offset_for(spectrum_name)
        buffer = self._get_lines_for(offset)
        spectrum = self._parse(buffer, spectrum_number)
        return spectrum

