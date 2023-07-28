import re
import io
import os
import logging
import itertools
import warnings

from typing import (
    Any, Callable, Collection,
    Dict, List, Mapping, Optional,
    Set, Tuple, Iterable, DefaultDict)

from pyteomics import proforma

from mzlib import annotation

from mzlib.analyte import FIRST_ANALYTE_KEY, FIRST_INTERPRETATION_KEY, Analyte, ProteinDescription
from mzlib.spectrum import Spectrum, SPECTRUM_NAME
from mzlib.attributes import AttributeManager, AttributeSet, Attributed

from .base import (
    DEFAULT_VERSION, FORMAT_VERSION_TERM, _PlainTextSpectralLibraryBackendBase,
    LIBRARY_NAME_TERM, AttributeSetTypes, SpectralLibraryBackendBase,
    SpectralLibraryWriterBase)
from .utils import try_cast, open_stream, CaseInsensitiveDict


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# TODO: Name: could be Compound or SpectrumName
leader_terms = {
    "Name": SPECTRUM_NAME,
    "NAME": SPECTRUM_NAME,
    "Compound": SPECTRUM_NAME,
    "COMPOUND": SPECTRUM_NAME,
}


def _generate_numpeaks_keys():
    w1 = "num"
    w2 = "peaks"
    seps = (" ", "")
    cases = (str.lower, str.title, str.upper)
    w1_cases = [c(w1) for c in cases]
    w2_cases = [c(w2) for c in cases]
    return {(w1c + sep + w2c) for (sep, w1c, w2c) in
            itertools.product(seps, w1_cases, w2_cases)}

NUM_PEAKS_KEYS = _generate_numpeaks_keys()

leader_terms_pattern = re.compile(r"(Name|NAME|Compound|COMPOUND)\s*:")
leader_terms_line_pattern = re.compile(r'(?:Name|NAME|Compound|COMPOUND)\s*:\s+(.+)')

STRIPPED_PEPTIDE_TERM = "MS:1000888|stripped peptide sequence"
PEPTIDE_MODIFICATION_TERM = "MS:1001471|peptide modification details"

PEAK_OBSERVATION_FREQ = "MS:1003279|observation frequency of peak"
PEAK_ATTRIB = "MS:1003254|peak attribute"


class AttributeHandler:
    keys: Collection[str]

    def __init__(self, keys: Collection[str]):
        if isinstance(keys, str):
            keys = (keys, )
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
            self.add_value(trans_key, try_cast(value), container)
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
                        self.add_group(k, map(try_cast, v), container)
            else:
                return False
        elif isinstance(trans_key, AttributeHandler):
            return trans_key(key, value, container)
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

    def __init__(self, chain: List[AttributeHandler] = None):
        if not chain:
            chain = []
        super().__init__(chain)
        self.mapping = CaseInsensitiveDict()
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



analyte_terms = CaseInsensitiveDict({
    "Charge": "MS:1000041|charge state",
    "precursor_charge": "MS:1000041|charge state",
    "precursorcharge": "MS:1000041|charge state",

    "MW": "MS:1000224|molecular mass",
    "total exact mass": "MS:1000224|molecular mass",
    "ExactMass": "MS:1000224|molecular mass",
    "exact_mass": "MS:1000224|molecular mass",
    "exact mass": "MS:1000224|molecular mass",
    "molecular formula": "MS:1000866|molecular formula",
    "Formula": "MS:1000866|molecular formula",
    "formula": "MS:1000866|molecular formula",
    "SMILES": "MS:1000868|SMILES formula",
    "InChIKey": "MS:1002894|InChIKey",
    # "Theo_mz_diff": "MS:1003209|monoisotopic m/z deviation",
    "Scan": {
        "Mods": PEPTIDE_MODIFICATION_TERM,
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
    "Mods": PEPTIDE_MODIFICATION_TERM,
    "Naa": "MS:1003043|number of residues",
    "PrecursorMonoisoMZ": "MS:1003208|experimental precursor monoisotopic m/z",
    "Mz_exact": "MS:1003208|experimental precursor monoisotopic m/z",
    "Mz_av": "MS:1003054|theoretical average m/z",
})


_HCD = ["MS:1000044|dissociation method", "MS:1000422|beam-type collision-induced dissociation"]

# TODO: qtof -> CAD
instrument_dispatch = CaseInsensitiveDict({
    "it": [["MS:1000044|dissociation method", "MS:1002472|trap-type collision-induced dissociation"]],
    "hcd": [_HCD],
    "QExactive": [_HCD],
    "Elite": [_HCD],
})


other_terms = CaseInsensitiveDict({
    "Parent": "MS:1003208|experimental precursor monoisotopic m/z",
    "ObservedPrecursorMZ": "MS:1003208|experimental precursor monoisotopic m/z",
    "PrecursorMZ": "MS:1003208|experimental precursor monoisotopic m/z",
    "PRECURSORMZ": "MS:1003208|experimental precursor monoisotopic m/z",
    "precursor": "MS:1003208|experimental precursor monoisotopic m/z",
    "precursor_mass": "MS:1003208|experimental precursor monoisotopic m/z",
    "precursormass": "MS:1003208|experimental precursor monoisotopic m/z",

    "Single": ["MS:1003065|spectrum aggregation type", "MS:1003066|singleton spectrum"],
    "Consensus": ["MS:1003065|spectrum aggregation type", "MS:1003067|consensus spectrum"],
    "Inst": instrument_dispatch,
    "Instrument_type": instrument_dispatch,
    "Spec": {"Consensus": [["MS:1003065|spectrum aggregation type", "MS:1003067|consensus spectrum"]]},
    "Scan": "MS:1003057|scan number",
    "Origfile": "MS:1003203|constituent spectrum file",
    "filename": "MS:1003203|constituent spectrum file",
    "file_name": "MS:1003203|constituent spectrum file",
    "Sample": "MS:1000002|sample name",
    "Filter": "MS:1000512|filter string",
    "FTResolution": "MS:1000028|detector resolution",
    "ms1PrecursorAb": "MS:1003085|previous MS1 scan precursor intensity",
    "Precursor1MaxAb": "MS:1003086|precursor apex intensity",
    "Purity": "MS:1009013|isolation window precursor purity",
    "Num peaks": "MS:1003059|number of peaks",
    "Num Peaks": "MS:1003059|number of peaks",
    "num_peaks": "MS:1003059|number of peaks",
    "numpeaks": "MS:1003059|number of peaks",
    "Run": "MS:1003203|constituent spectrum file",
    "Splash": "MS:1002599|splash key",
})


@FunctionAttributeHandler.wraps("num_unassigned_peaks")
def unassigned_peaks_handler(key: str, value: str, container: Attributed) -> bool:
    is_top_20 = False
    if isinstance(value, str):
        if "/" in value:
            if "/20" in value:
                is_top_20 = True
            value = value.split("/")[0]
        value = int(value)
    if is_top_20:
        container.add_attribute("MS:1003290|number of unassigned peaks among top 20 peaks", value)
    else:
        container.add_attribute("MS:1003288|number of unassigned peaks", value)
    return True


interpretation_terms = CaseInsensitiveDict({
    "Unassigned_all_20ppm": "MS:1003079|total unassigned intensity fraction",
    "Unassign_all": "MS:1003079|total unassigned intensity fraction",

    "top_20_num_unassigned_peaks_20ppm": "MS:1003290|number of unassigned peaks among top 20 peaks",
    "num_unassigned_peaks_20ppm": unassigned_peaks_handler,
    "num_unassigned_peaks": unassigned_peaks_handler,

    "max_unassigned_ab_20ppm": "MS:1003289|intensity of highest unassigned peak",
    "max_unassigned_ab": "MS:1003289|intensity of highest unassigned peak",

    "Unassigned_20ppm": "MS:1003080|top 20 peak unassigned intensity fraction",
    "Unassigned": "MS:1003080|top 20 peak unassigned intensity fraction",
})


interpretation_member_terms = CaseInsensitiveDict({
    "Q-value": "MS:1002354|PSM-level q-value",
})


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
    "rat": [
        ["MS:1001467|taxonomy: NCBI TaxID", "NCBITaxon:10116|Rattus norvegicus"],
        ["MS:1001469|taxonomy: scientific name",
         "Rattus norvegicus"],
        ["MS:1001468|taxonomy: common name", "rat"]
    ]
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

TERMINAL_MODIFICATIONS = {
    "Acetyl",
    "TMT6plex"
}


# TODO: ppm is unsigned, add mass calculation to determine true mass accuracy

annotation_pattern = re.compile(r"""^
(?:(?:(?P<series>[abyxcz]\.?)(?P<ordinal>\d+))|
   (:?(?P<series_internal>(?:Int/[ARNDCEQGHKMFPSTWYVILJarndceqghkmfpstwyvilj]+)|(?:[mM](?P<internal_start>\d+):(?P<internal_end>\d+))))|
   (?P<precursor>p)|
   (:?I(?P<immonium>[ARNDCEQGHKMFPSTWYVIL])(?:(?P<immonium_modification>CAM)|[A-Z])?)|
   (?P<reporter>(?P<reporter_label>TMT|ITRAQ|iTRAQ)(?P<reporter_mass>\d+[NC]?))|
   (?:_(?P<external_ion>[^\s,/]+))|
   (?P<unannotated>\?)
)
(?P<neutral_losses>((?:[+-]\d*[A-Z][A-Za-z0-9]*)|(?:[+-]iTRAQ|TMT)|(?:[+-]\d+?(?!i)))+)?
(?:(?:\[M(?P<adducts>(:?[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?
(?:\^(?P<charge>[+-]?\d+))?
(?:(?P<isotope>[+-]?\d*)i)?)+
(?:/(?:Dev=)?(?P<mass_error>[+-]?\d+(?:\.\d+))(?P<mass_error_unit>ppm)?)?
""", re.X)


class MSPAnnotationStringParser(annotation.AnnotationStringParser):
    def _dispatch_internal_peptide_fragment(self, data: Dict[str, Any], adducts: List, charge: int, isotope: int, neutral_losses: List,
                                            analyte_reference: Any, mass_error: Any, **kwargs):
        if data['internal_start']:
            # The hybrid pattern detected an m<start>:<end> type string, not an Int/seq string
            return super(MSPAnnotationStringParser, self)._dispatch_internal_peptide_fragment(
                data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, **kwargs)

        spectrum = kwargs.get("spectrum")
        if spectrum is None:
            raise ValueError("Cannot infer sequence coordinates from MSP internal fragmentation notation without"
                             " a reference to the spectrum, please pass spectrum as a keyword argument")
        sequence = self._get_peptide_sequence_for_analyte(
            spectrum, analyte_reference)
        subseq = data['series_internal']
        if subseq.startswith("Int/"):
            subseq = subseq[4:]
        subseq = subseq.upper()
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

    def _coerce_isotope(self, data):
        value = data.get('isotope')
        if value is not None:
            if value == '':
                data['isotope'] = '1'
        return super()._coerce_isotope(data)

    def _coerce_neutral_losses(self, data: Dict[str, str]) -> List:
        return super()._coerce_neutral_losses(data)

    def _coerce_analyte_reference(self, data: Dict[str, str]) -> str:
        return None

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


def null_handler(key: str, value: str, container: Attributed) -> bool:
    return True


msp_spectrum_attribute_handler = DispatchingAttributeHandler()
msp_analyte_attribute_handler = DispatchingAttributeHandler()


msp_spectrum_attribute_handler.add(FunctionAttributeHandler("Nprot", null_handler))
msp_spectrum_attribute_handler.add(FunctionAttributeHandler("Peptype", null_handler))


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("Spectrum_type")
def ms_level_handler(key: str, value: str, container: Attributed) -> bool:
    attr_name = "MS:1000511|ms level"
    if value is None:
        return False
    if isinstance(value, str):
        if value.lower().startswith("ms"):
            value = int(value[2:])
    container.add_attribute(attr_name, value)
    return True


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("ionmode", "ion_mode", "ionization mode", "IONMODE", "Ion_mode")
def polarity_handler(key: str, value: str, container: Attributed) -> bool:
    polarity_term = "MS:1000465|scan polarity"
    positive = "MS:1000130|positive scan"
    negative = "MS:1000129|negative scan"

    if isinstance(value, str):
        value = value.lower()
        if value == 'positive':
            value = positive
        elif value == 'negative':
            value = negative
        else:
            raise ValueError(f"Can't infer scan polarity from {value}")
        container.add_attribute(polarity_term, value)
        return True
    return False


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
@FunctionAttributeHandler.wraps("Collision_energy", "CE", "colenergy", "collisionenergy", "ionization energy")
def collision_energy_handler(key: str, value: str, container: Attributed) -> bool:
    if isinstance(value, str):
        if "NCE" in value:
            return normalized_collision_energy_handler(key, value, container)
        match = re.match(r"([\d\.]+)", value)
        if match is not None:
            value = try_cast(match.group(1))
    if value is not None:
        group_identifier = container.get_next_group_identifier()
        container.add_attribute(
            "MS:1000045|collision energy", value, group_identifier)
        container.add_attribute(
            "UO:0000000|unit", "UO:0000266|electronvolt", group_identifier)
        return True
    return False


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("NCE", "nce")
def normalized_collision_energy_handler(key: str, value: str, container: Attributed) -> bool:
    if isinstance(value, str):
        match = re.match(r"([\d\.]+)", value)
        if match is not None:
            value = try_cast(match.group(1))
    if value is not None:
        group_identifier = container.get_next_group_identifier()
        container.add_attribute(
            "MS:1000138|normalized collision energy", value, group_identifier)
        container.add_attribute(
            "UO:0000000|unit", "UO:0000187|percent", group_identifier)
        return True
    return False


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("RT", "rettime", "retentiontime", "rtinseconds")
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
            if float(match.group(1)) > 250 or key.lower() == 'rtinseconds':
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
            "MS:1001975|delta m/z", abs(value), group_identifier)
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
def protein_handler(key, value, container: Attributed):
    if value is None:
        return False
    key = "MS:1000885|protein accession"
    match = re.match(r"\(pre=(.),post=(.)\)", value)
    group_identifier = None
    if match is not None:
        value = value[:match.start()]
        group_identifier = container.get_next_group_identifier()
        container.add_attribute("MS:1001112|n-terminal flanking residue",
                                match.group(1), group_identifier=group_identifier)
        container.add_attribute("MS:1001113|c-terminal flanking residue",
                                match.group(2), group_identifier=group_identifier)
    container.add_attribute(key, re.sub(r"\(pre=(.),post=(.)\)", '', value.strip('"').strip("'")),
                            group_identifier=group_identifier)
    return True


@msp_spectrum_attribute_handler.add
@FunctionAttributeHandler.wraps("BasePeak")
def base_peak_handler(key, value, container: Attributed):
    if value is None:
        return False
    value = float(value)
    group_id = container.get_next_group_identifier()
    container.add_attribute("MS:1000505|base peak intensity", value, group_id)
    container.add_attribute("UO:0000000|unit", "MS:1000131|number of detector counts", group_id)
    return True


class _UnknownTermTracker:
    counts: DefaultDict

    def add(self, key: str, value: Optional[str]=None):
        raise NotImplementedError()

    def items(self):
        return self.counts.items()


class UnknownKeyValueTracker(_UnknownTermTracker):
    def __init__(self) -> None:
        self.counts = DefaultDict(lambda: DefaultDict(int))

    def add(self, key: str, value: Optional[str]=None):
        self.counts[key][value] += 1


class UnknownKeyTracker(_UnknownTermTracker):
    def __init__(self) -> None:
        self.counts = DefaultDict(int)

    def add(self, key: str, value: Optional[str] = None):
        self.counts[key] += 1


protein_attributes_to_group = [
    "MS:1003048|number of enzymatic termini",
    "MS:1001045|cleavage agent name",
    "MS:1001112|n-terminal flanking residue",
    "MS:1001113|c-terminal flanking residue",
    "MS:1003044|number of missed cleavages",
    "MS:1000885|protein accession",
]


class MSPSpectralLibrary(_PlainTextSpectralLibraryBackendBase):
    file_format = "msp"
    format_name = "msp"

    modification_parser: ModificationParser
    unknown_attributes: _UnknownTermTracker

    def __init__(self, filename, index_type=None, read_metadata=True, create_index: bool=True):
        super().__init__(filename, index_type, read_metadata, create_index=create_index)
        self.modification_parser = ModificationParser()
        self.unknown_attributes = UnknownKeyTracker()

    @classmethod
    def guess_from_header(cls, filename: str) -> bool:
        with open_stream(filename, 'r') as stream:
            first_line = stream.readline()
            if leader_terms_pattern.match(first_line):
                return True
        return False

    def read_header(self) -> bool:
        filename = self.filename
        file_like_object = isinstance(filename, io.IOBase)

        stream = open_stream(filename)
        match, offset = self._parse_header_from_stream(stream)
        if file_like_object:
            if stream.seekable():
                stream.seek(0)
            stream.detach()
        else:
            stream.close()
        return match

    def _parse_header_from_stream(self, stream: io.IOBase) -> Tuple[bool, int]:
        first_line = stream.readline()
        attributes = AttributeManager()
        attributes.add_attribute(FORMAT_VERSION_TERM, DEFAULT_VERSION)
        if isinstance(self.filename, (str, os.PathLike)):
            attributes.add_attribute(LIBRARY_NAME_TERM, self.filename.rsplit('.msp', 1)[0].split(os.sep)[-1])
        elif hasattr(stream, 'name'):
            attributes.add_attribute(LIBRARY_NAME_TERM, stream.name.rsplit('.msp', 1)[0].split(os.sep)[-1])
        self.attributes.clear()
        self.attributes._from_iterable(attributes)
        if leader_terms_pattern.match(first_line):
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

        begin_loc = None
        file_like_object = False
        #### Determine the filesize
        try:
            file_size = os.path.getsize(filename)
        except TypeError:
            if isinstance(filename, io.IOBase):
                file_like_object = True
                if filename.seekable():
                    begin_loc = filename.tell()
                    filename.seek(0, os.SEEK_END)
                    file_size = filename.tell()
                    filename.seek(0)
        infile = open_stream(filename, 'r')

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
            # TODO: Name: could be Compound or SpectrumName
            if state == 'header':
                if leader_terms_pattern.match(line):
                    state = 'body'
                    spectrum_file_offset = line_beginning_file_offset
                else:
                    continue
            if state == 'body':
                if len(line) == 0:
                    continue
                if leader_terms_pattern.match(line):
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
                    spectrum_name = leader_terms_line_pattern.match(line).group(1)

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

        if not file_like_object:
            infile.close()
        else:
            infile.detach()
            if begin_loc is not None:
                filename.seek(begin_loc)
        return n_spectra

    def _buffer_from_stream(self, infile: Iterable[str]) -> List[str]:
        state = 'body'
        spectrum_buffer = []

        file_offset = 0
        line_beginning_file_offset = 0
        for line in infile:
            line_beginning_file_offset = file_offset
            file_offset += len(line)
            line = line.rstrip()
            if state == 'body':
                if len(line) == 0:
                    continue
                if leader_terms_pattern.match(line):
                    if len(spectrum_buffer) > 0:
                        return spectrum_buffer
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
                    attributes[key] = value
                elif line.count("=") > 0:
                    key, value = re.split(r"=\s*", line, 1)
                    attributes[key] = value
                elif line.count("\t") > 0:
                    warnings.warn(f"Line {line!r} looks like a peak annotation?")
                    in_header = False
                else:
                    key = line
                    value = None
                    attributes[key] = value


                #### If the key is "Num peaks" then we're done with the header and peak list follows
                if key in NUM_PEAKS_KEYS:
                    in_header = False

                #### The "Comment" key requires special parsing
                if key == "Comment" or key == "Comments":

                    #### Remove it from attributes
                    del attributes[key]
                    self._parse_comment(value, attributes)

            #### Else in the peaks section. Parse the peaks.
            else:
                #### Split into the expected three values
                values = re.split(r'\s+', line, maxsplit=2)
                # Sometimes MSP files have multiple peaks on the same line, delimited by ';',
                # so we must potentially parse more than one peak per line.
                if values[1].endswith(";"):
                    buffered_peaks = []
                    buffered_peaks.append(values[:2])
                    buffered_peaks[0][1] = buffered_peaks[0][1].strip(";")
                    rest = values[0]
                    for block in re.split(r";\s?", rest):
                        if block:
                            buffered_peaks.append(re.split(r'\s+', block, maxsplit=2))

                else:
                    buffered_peaks = [values]

                for values in buffered_peaks:
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
        has_observation_freq = (False, 0)
        for i, peak in enumerate(spectrum.peak_list):
            try:
                parsed_interpretation = self._parse_annotation(peak[2], spectrum=spectrum)
            except ValueError as err:
                message = str(err)
                raise ValueError(
                    f"An error occurred while parsing the peak annotation for peak {i}: {message}") from err

            if parsed_interpretation and isinstance(parsed_interpretation[0], annotation.InvalidAnnotation) and peak[3]:
                logger.debug("Failed to parse interpretation string %r with assumed %d aggregation fields, trying to parse the combined string",
                             peak[2], len(peak[3]))

                peak[2] = " ".join([peak[2]] + peak[3])
                peak[3] = []
                try:
                    parsed_interpretation = self._parse_annotation(peak[2], spectrum=spectrum)
                except ValueError as err:
                    message = str(err)
                    raise ValueError(
                        f"An error occurred while parsing the peak annotation for peak {i}: {message}") from err

            peak[2] = parsed_interpretation
            for i, agg in enumerate(peak[3]):
                # NOTE: Making the potentially unsafe assumption that the only fractional
                # aggregation statistic is the peak frequency. At this time, we haven't formalized
                # any other aggregation metric that I also have an MSP example for.
                if '/' in agg:
                    try:
                        agg = _parse_fraction(agg)
                        peak[3][i] = agg
                        has_observation_freq = (True, i)
                    except ValueError:
                        pass
                    except Exception as err:
                        raise ValueError(f"An error occurred while parsing the peak aggregation for peak {i}: {str(err)}") from err

        aggregation_metrics = []
        if has_observation_freq[0]:
            aggregation_metrics.append((PEAK_ATTRIB, PEAK_OBSERVATION_FREQ))

        if aggregation_metrics:
            spectrum.add_attribute_group(aggregation_metrics)

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
            if n_quotes % 2 == 0:
                fixed_comment_items.append(new_item)
                new_item = ""

        #### Try to split each item on the first = character and store in attributes
        for item in fixed_comment_items:
            #### If there is an = then split on the first one and store key and value
            if item.count("=") > 0:
                comment_key, comment_value = item.split("=", 1)
                cleaned_key = comment_key.strip('"')
                if len(cleaned_key) != len(comment_key):
                    comment_value = comment_value.strip('"')
                attributes[cleaned_key] = try_cast(comment_value)

            #### Otherwise just store the key with a null value
            else:
                attributes[item] = None

    def _make_attribute_handlers(self):
        other_manager = MappingAttributeHandler(other_terms)
        analyte_manager = MappingAttributeHandler(analyte_terms)
        interpretation_manager = MappingAttributeHandler(interpretation_terms)
        interpretation_member_manager = MappingAttributeHandler(interpretation_member_terms)
        return (other_manager,
                analyte_manager,
                interpretation_manager,
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
                break
        else:
            spectrum.add_attribute(
                "ERROR", f"Required term {leader_terms[term]} is missing")

        #### Translate the rest of the known attributes and collect unknown ones
        unknown_terms = []

        (other_manager, analyte_manager, interpretation_manager, interpretation_member_manager,
         spectrum_func_handler, analyte_func_handler) = self._make_attribute_handlers()
        for attribute in attributes:

            #### Skip a leader term that we already processed
            if attribute in leader_terms:
                continue
            if not attribute:
                continue
            if attribute in other_manager:
                if not other_manager(attribute, attributes[attribute], spectrum):
                    unknown_terms.append(attribute)

            elif attribute in analyte_manager:
                if not analyte_manager(attribute, attributes[attribute], analyte):
                    unknown_terms.append(attribute)

            elif attribute in interpretation_manager:
                if not interpretation_manager(attribute, attributes[attribute], interpretation):
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
                            analyte.add_attribute(
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
                if name:
                    match = re.match(r"([ARNDCEQGHKMFPSTWYVIL]+)/(\d+)", name)
                    if match:
                        analyte.add_attribute(
                            STRIPPED_PEPTIDE_TERM, match.group(1))
                        analyte.add_attribute(
                            "MS:1000041|charge state", try_cast(match.group(2)))

        #### Handle the uninterpretable terms
        for attribute in unknown_terms:
            self.unknown_attributes.add(attribute, attributes[attribute])
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

    def _deduplicate_attr(self, container: Attributed, attribute_name: str) -> int:
        try:
            attrs = container.get_attribute(attribute_name, raw=True)
        except KeyError:
            return 0
        if isinstance(attrs, list):
            # Sometimes there are multiple entries that map to molecular mass, and one is the nominal
            # mass. Prefer listing the most precise mass first. Apply the same logic to all attributes.
            attrs.sort(key=lambda x: len(str(x.value)), reverse=True)
            container.remove_attribute(attrs[0].key)
            seen_values = set()
            for attr in attrs:
                key = (attr.value, attr.group_id)
                if key in seen_values:
                    continue
                seen_values.add(key)
                container.add_attribute(attr.key, attr.value, attr.group_id)
            return len(attrs) - len(seen_values)
        return 0

    def _complete_analyte(self, analyte: Analyte):
        if analyte.has_attribute(STRIPPED_PEPTIDE_TERM):
            peptide = proforma.ProForma.parse(analyte.get_attribute(STRIPPED_PEPTIDE_TERM))
            if analyte.has_attribute(PEPTIDE_MODIFICATION_TERM):
                modification_details = analyte.get_attribute(PEPTIDE_MODIFICATION_TERM)
                mods = self.modification_parser(modification_details)
                for position, residue, mod in mods:
                    if position == 0 and mod in TERMINAL_MODIFICATIONS:
                        peptide.n_term = [proforma.GenericModification(mod)]
                    else:
                        seqpos = list(peptide.sequence[position])
                        if not seqpos[1]:
                            seqpos[1] = [proforma.GenericModification(mod)]
                        peptide.sequence[position] = tuple(seqpos)
                        assert seqpos[0] == residue
            analyte.add_attribute("MS:1003169|proforma peptidoform sequence", str(peptide))
            if analyte.has_attribute(PEPTIDE_MODIFICATION_TERM):
                analyte.remove_attribute(PEPTIDE_MODIFICATION_TERM)
            analyte.add_attribute("MS:1001117|theoretical mass", peptide.mass)

        self._pack_protein_description(analyte)

        self._deduplicate_attr(analyte, "MS:1000224|molecular mass")
        self._deduplicate_attr(analyte, "MS:1000866|molecular formula")

    def _pack_protein_description(self, analyte: Analyte):
        table = {}
        max_len = 0

        # Collect terms that describe a protein
        for term in protein_attributes_to_group:
            if not analyte.has_attribute(term):
                table[term] = []
                continue
            values = analyte.get_attribute(term, raw=True)
            if not isinstance(values, list):
                values = [values]
            max_len = max(max_len, len(values))
            table[term] = values

        # Ensure that all arrays are the same length
        for k, v in table.items():
            if len(v) < max_len:
                v.extend([None] * (max_len - len(v)))

        # Group together terms and remove the previous entries
        groups = []
        for i in range(max_len):
            group = []
            for k, v in table.items():
                if v[i] is None:
                    continue
                group.append((k, v[i].value))
                analyte.remove_attribute(v[i].key, v[i].group_id)
            groups.append(group)

        # Now add back the groups
        for group in groups:
            analyte.add_attribute_group(group)

    def get_spectrum(self, spectrum_number: int=None, spectrum_name: str=None) -> Spectrum:
        # keep the two branches separate for the possibility that this is not possible with all
        # index schemes.
        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError("Provide only one of spectrum_number or spectrum_name")
            index_record = self.index.record_for(spectrum_number)
            offset = index_record.offset
        elif spectrum_name is not None:
            index_record = self.index.record_for(spectrum_name)
            spectrum_number = index_record.number
            offset = index_record.offset
        buffer = self._get_lines_for(offset)
        spectrum = self._parse(buffer, index_record.index)
        return spectrum

    def summarize_parsing_errors(self) -> Dict:
        errors = super().summarize_parsing_errors()
        errors['unknown attributes'] = DefaultDict(dict)
        for k, v in self.unknown_attributes.items():
            errors['unknown attributes'][k] = v
        return errors


def _parse_fraction(x: str) -> float:
    a, b = x.split("/")
    return int(a) / int(b)


class MSPSpectralLibraryWriter(SpectralLibraryWriterBase):
    file_format = "msp"

    analyte_keys = {
        "MS:1000866|molecular formula": "Formula",
        "MS:1003044|number of missed cleavages": "MC",
        "MS:1001471|peptide modification details": "Mods",
        "MS:1003043|number of residues": "Naa",
        "MS:1003208|experimental precursor monoisotopic m/z": "PrecursorMonoisoMZ",
        "MS:1003054|theoretical average m/z": "Mz_av",
        "MS:1003169|proforma peptidoform sequence": "ProForma",
        "MS:1000888|stripped peptide sequence": "Peptide",
        "MS:1000041|charge state": "Charge",
    }

    for species_name, keys in species_map.items():
        analyte_keys[tuple(keys[0])] = ("Organism", species_name)

    modification_map = {v: k for k, v in MODIFICATION_NAME_MAP.items()}

    spectrum_keys = {
        "MS:1000041|charge state": "Charge",
        ("MS:1003065|spectrum aggregation type", "MS:1003066|singleton spectrum"): "Single",
        ("MS:1003065|spectrum aggregation type", "MS:1003067|consensus spectrum"): "Consensus",
        "MS:1003057|scan number": "Scan",
        "MS:1003203|constituent spectrum file": "Origfile",
        "MS:1000002|sample name": "Sample",
        "MS:1000512|filter string": "Filter",
        "MS:1003086|precursor apex intensity": "Precursor1MaxAb",
        "MS:1009013|isolation window precursor purity": "Purity",
        "MS:1000505|base peak intensity": "BasePeak",
        "MS:1002599|splash key": "Splash",
        "MS:1003289|intensity of highest unassigned peak": "max_unassigned_ab",
        "MS:1003080|top 20 peak unassigned intensity fraction": "Unassigned",
        "MS:1003079|total unassigned intensity fraction": "Unassign_all",
        "MS:1003290|number of unassigned peaks among top 20 peaks": "top_20_num_unassigned_peaks",
        ("MS:1000044|dissociation method", "MS:1000422|beam-type collision-induced dissociation"): "HCD",
        ("MS:1000044|dissociation method", "MS:1002472|trap-type collision-induced dissociation"): "CID",
    }

    # TODO: add these
    interpretation_keys = {
        "MS:1002354|PSM-level q-value": "Q-value",
    }

    def __init__(self, filename, **kwargs):
        super(MSPSpectralLibraryWriter, self).__init__(filename)
        self._coerce_handle(self.filename)

    def write_header(self, library: SpectralLibraryBackendBase):
        pass

    def write_attribute_set(self, attribute_set: AttributeSet, attribute_set_type: AttributeSetTypes):
        pass

    def _format_value(self, value):
        if isinstance(value, str):
            if not (value.startswith('"') and value.endswith('"')):
                value = f"\"{value}\""
        return str(value)

    def _proforma_to_mods(self, proforma_seq: str) -> str:
        parsed = proforma.ProForma.parse(proforma_seq)
        mods = [(i, tok) for i, tok in enumerate(parsed) if tok[1]]
        if mods:
            tokens = []
            for i, mod_site in mods:
                tokens.append(str(i))
                tokens.append(mod_site[0])
                mod = mod_site[1][0]
                mod_name = self.modification_map.get(mod.name, mod.name)
                tokens.append(mod_name)
            return f'Mods={len(mods)}({",".join(tokens)})'
        else:
            return 'Mods=0'

    def _protein_to_comments(self, analyte: Analyte) -> List[str]:
        acc = []
        protein: ProteinDescription
        for protein in analyte.proteins:
            accession = None
            pre = None
            post = None
            if protein.accession:
                accession = protein.accession
            if protein.flanking_n_terminal_residue:
                pre = protein.flanking_n_terminal_residue
            if protein.flanking_c_terminal_residue:
                post = protein.flanking_c_terminal_residue
            if accession:
                token = accession
                if token.startswith('"'):
                    token = token.strip('"')
                if pre or post:
                    token += f"(pre={pre or '-'},post={post or '-'})"
                acc.append(f"Protein={self._format_value(token)}")
                if protein.number_of_enzymatic_termini == 2:
                    if protein.cleavage_agent == "MS:1001251|Trypsin":
                        acc.append("Pep=Tryptic")
                elif protein.number_of_enzymatic_termini == 1:
                    if protein.cleavage_agent == "MS:1001251|Trypsin":
                        acc.append("Pep=SemiTryptic")
                if protein.missed_cleavages is not None:
                    acc.append(f"MC={self._format_value(protein.missed_cleavages)}")
                break
        return acc

    def _build_comments(self, spectrum: Spectrum, attribute_container: Attributed,
                        rule_map: Dict) -> List[Tuple[str, str]]:
        accumulator = []

        for attr_name, msp_name in rule_map.items():
            if isinstance(attr_name, tuple):
                if attribute_container.has_attribute(attr_name[0]):
                    value = attribute_container.get_attribute(attr_name[0])
                    if value == attr_name[1]:
                        if isinstance(msp_name, str):
                            accumulator.append(msp_name)
                        elif isinstance(msp_name, (list, tuple)) and len(msp_name) == 2:
                            accumulator.append('='.join(msp_name))
                        else:
                            raise TypeError(f"Can't infer conversion for {msp_name} given {attr_name}")
            elif attribute_container.has_attribute(attr_name):
                value = attribute_container.get_attribute(attr_name)
                if isinstance(value, list):
                    logger.warn(
                        "Spectrum %r contains multiple values for %r, only the first will be saved",
                        spectrum.name, attr_name
                    )
                    accumulator.append(f"{msp_name}={self._format_value(value[0])}")
                else:
                    accumulator.append(f"{msp_name}={self._format_value(value)}")
        return accumulator

    def build_spectrum_comments(self, spectrum: Spectrum) -> List[Tuple[str, str]]:
        accumulator = self._build_comments(spectrum, spectrum, self.spectrum_keys)
        if spectrum.analytes:
            analyte = spectrum.get_analyte('1')
            accumulator += self._build_comments(spectrum, analyte, self.analyte_keys)
            if analyte.peptide:
                accumulator.append(self._proforma_to_mods(analyte.peptide))
            accumulator += self._protein_to_comments(analyte)
        if spectrum.interpretations:
            interp = spectrum.get_interpretation('1')
            accumulator += self._build_comments(
                spectrum, interp, self.interpretation_keys)
        return accumulator

    def write_spectrum(self, spectrum: Spectrum):
        if len(spectrum.analytes) > 1:
            logger.warning(
                "Spectrum %r contains multiple analytes, MSP will only contain the first",
                spectrum.name
            )
        analyte = spectrum.get_analyte('1')
        self.handle.write(f"Name: {spectrum.name}\n")
        self.handle.write(f"MW: {analyte.mass}\n")
        self.handle.write(f"Comment: {' '.join(self.build_spectrum_comments(spectrum))}\n")
        self._write_peaks(spectrum)
        self.handle.write("\n")

    def _format_annotation(self, annot: annotation.IonAnnotationBase):
        parts = []
        if isinstance(annot, annotation.PeptideFragmentIonAnnotation):
            parts.append(f"{annot.series}{annot.position}")
        elif isinstance(annot, annotation.ImmoniumIonAnnotation):
            parts.append(f"I{annot.amino_acid}{annot.modification if annot.modification else ''}")
        elif isinstance(annot, annotation.ReporterIonAnnotation):
            parts.append(annot.reporter_label)
        if not parts:
            return "?"
        if annot.neutral_losses:
            f = annotation.combine_formula(annot.neutral_losses)
            if f[0] not in ('-', '+'):
                f = '+' + f
            parts.append(f)
        if annot.adducts:
            parts.append(f"[{annotation.combine_formula(annot.adducts)}]")
        if annot.charge > 1:
            parts.append(f"^{annot.charge}")
        if annot.isotope:
            if annot.isotope > 0:
                parts.append(f'+{annot.isotope}i')
            else:
                parts.append(f'{annot.isotope}i')
        if annot.mass_error:
            parts.append(f'/{annot.mass_error}')
        return ''.join(parts)

    def _write_peaks(self, spectrum: Spectrum):
        self.handle.write(f"Num peaks: {len(spectrum.peak_list)}\n")
        for peak in spectrum.peak_list:
            annot = self._format_annotation(peak[2][0]) if peak[2] else '?'
            self.handle.write(f"{peak[0]:0.4f}\t{peak[1]:0.4f}\t\"{annot}\"\n")

    def close(self):
        self.handle.close()
