import re
from sys import intern
from typing import Any, List, Pattern, Dict, Tuple, Union


JSONDict = Dict[str, Union[List, Dict, int, float, str, bool, None]]

annotation_pattern = re.compile(r"""
^(?P<is_auxiliary>&)?
   (?:(?P<analyte_reference>\d+)@)?
   (?:(?:(?P<series>[axbycz]\.?)(?P<ordinal>\d+))|
   (?P<series_internal>[m](?P<internal_start>\d+):(?P<internal_end>\d+))|
   (?P<precursor>p)|
   (:?I(?P<immonium>[ARNDCEQGHKMFPSTWYVIL])(?:\[(?P<immonium_modification>(?:[^\]]+))\])?)|
   (?P<reporter>r(?:
    (?:\[
        (?P<reporter_label>[^\]]+)
    \])
   ))|
   (?:f\{(?P<formula>[A-Za-z0-9]+)\})|
   (?:_\{
    (?P<external_ion>[^\{\}\s,/]+)
    \})|
   (?:s\{(?P<smiles>[^\}]+)\})|
   (?:(?P<unannotated>\?)(?P<unannotated_label>\d+)?)
)
(?P<neutral_losses>(?:[+-]\d*
    (?:(?:[A-Z][A-Za-z0-9]*)|
        (?:\[
            (?:
                (?:[A-Za-z0-9:\.]+)
            )
            \])
    )
)+)?
(?:(?P<isotope>[+-]\d*)i)?
(?:\^(?P<charge>[+-]?\d+))?
(?:\[(?P<adducts>M(:?[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?
(?:/(?P<mass_error>[+-]?\d+(?:\.\d+)?)(?P<mass_error_unit>ppm)?)?
(?:\*(?P<confidence>\d*(?:\.\d+)?))?
""", re.X)

# At the time of first writing, this pattern could be translated into the equivalent
# ECMAScript compliant regex:
# ^(?:(?<analyte_reference>[^@\s]+)@)?(?:(?:(?<series>[axbycz]\.?)(?<ordinal>\d+))|(?<series_internal>[m](?<internal_start>\d+):(?<internal_end>\d+))|(?<precursor>p)|(:?I(?<immonium>[ARNDCEQGHKMFPSTWYVIL])(?:\[(?<immonium_modification>(?:[^\]]+))\])?)|(?<reporter>r(?:(?:\[(?<reporter_label>[^\]]+)\])))|(?:f\{(?<formula>[A-Za-z0-9]+)\})|(?:_(?<external_ion>[^\s,/]+)))(?<neutral_losses>(?:[+-]\d*(?:(?:[A-Z][A-Za-z0-9]*)|(?:\[(?:(?:[A-Za-z0-9:\.]+))\])))+)?(?:(?<isotope>[+-]\d*)i)?(?:\^(?<charge>[+-]?\d+))?(?:\[(?<adducts>M(:?[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?(?:/(?<mass_error>[+-]?\d+(?:\.\d+)?)(?<mass_error_unit>ppm)?)?(?:\*(?<confidence>\d*(?:\.\d+)?))?
# Line breaks not introduced to preserve syntactic correctness.

def _sre_to_ecma(pattern):
    # Assumes that expected whitespace matches are denoted with \s
    return pattern.replace("?P<", "?<").replace("\n", '').replace(" ", "")


def tokenize_signed_symbol_list(string: str) -> List[str]:
    '''Parse a string containing signed lists of symbols into
    a signed token list

    Parameters
    ----------
    string: str
        The symbol sequence to parse

    Returns
    -------
    list
    '''
    if not string:
        return []
    tokens = re.split(r"\s*(\+|-)\s*", string)
    if not tokens:
        return []

    result = []
    sign = '+'
    symbol = None
    for token in tokens:
        if not token.strip():
            continue
        if token in ('+', '-'):
            if symbol is not None:
                result.append(sign + symbol if sign == '-' else symbol)
            sign = token
            symbol = None
        else:
            symbol = token
    if symbol is not None:
        result.append(sign + symbol if sign == '-' else symbol)
    return result


def combine_formula(tokens: List[str], leading_sign: bool=False):
    if not tokens:
        return ''
    if not tokens[0].startswith("-") and leading_sign:
        out = ['+' + tokens[0]]
    else:
        out = [tokens[0]]
    for token in tokens[1:]:
        if token.startswith("-") or token.startswith("+"):
            out.append(token)
        else:
            out.append('+' + token)
    return ''.join(out)


class MassError(object):
    _DEFAULT_UNIT = "Da"

    unit: str
    mass_error: float

    def __init__(self, mass_error, unit=None):
        if unit is None:
            unit = self._DEFAULT_UNIT
        self.mass_error = float(mass_error)
        self.unit = unit

    def serialize(self) -> str:
        unit = self.unit
        if unit == self._DEFAULT_UNIT:
            unit = ''
        return f"{self.mass_error}{unit}"

    def __eq__(self, other):
        return self.mass_error == other.mass_error and self.unit == other.unit

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return self.serialize()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.mass_error}, {self.unit})"

    def to_json(self) -> JSONDict:
        return {
            "value": self.mass_error,
            "unit": self.unit
        }


class SeriesLabelSubclassRegisteringMeta(type):
    def __new__(mcls, name, bases, attrs):
        label = attrs.get("series_label")
        if label and isinstance(label, str):
            label = intern(label)
        override_label = attrs.get('override_label', False)
        cls = super(SeriesLabelSubclassRegisteringMeta, mcls).__new__(mcls, name, bases, attrs)
        if not hasattr(cls, '_label_registry'):
            registry = cls._label_registry = {}
        else:
            registry = cls._label_registry
        if label:
            if label not in registry or override_label:
                registry[label] = cls
        return cls


class IonAnnotationBase(object, metaclass=SeriesLabelSubclassRegisteringMeta):
    __slots__ = ("series", "neutral_losses", "isotope", "adducts", "charge", "analyte_reference",
                 "mass_error", "confidence", "rest", "is_auxiliary")

    series_label = None
    _molecule_description_fields = {}

    series: str
    neutral_losses: List
    isotope: int
    adducts: List
    charge: int
    analyte_reference: str
    mass_error: MassError
    confidence: float
    rest: Any
    is_auxiliary: bool

    def __init__(self, series, neutral_losses=None, isotope=None, adducts=None, charge=None,
                 analyte_reference=None, mass_error=None, confidence=None, rest=None,
                 is_auxiliary=None):
        if isotope is None:
            isotope = 0
        if charge is None:
            charge = 1

        self.series = series
        self.neutral_losses = neutral_losses or []
        self.isotope = isotope
        self.adducts = adducts or []
        self.charge = charge
        self.analyte_reference = analyte_reference
        self.mass_error = mass_error
        self.confidence = confidence
        self.rest = rest
        self.is_auxiliary = is_auxiliary

    @property
    def adduct(self) -> List:
        return self.adducts

    @adduct.setter
    def adduct(self, value):
        self.adducts = value

    @property
    def neutral_loss(self) -> List:
        return self.neutral_losses

    @neutral_loss.setter
    def neutral_loss(self, value):
        self.neutral_losses = value

    def __hash__(self):
        return hash(self.serialize())

    def __eq__(self, other):
        return self.serialize() == str(other)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return self.serialize()

    def _format_ion(self) -> str:
        raise NotImplementedError()

    def serialize(self) -> str:
        parts = []
        if self.analyte_reference is not None:
            parts.append(f"{self.analyte_reference}@")
        parts.append(self._format_ion())
        if self.neutral_losses:
            parts.append(combine_formula(self.neutral_losses, leading_sign=True))
        if self.isotope != 0:
            sign = "+" if self.isotope > 0 else "-"
            isotope = abs(self.isotope)
            if isotope == 1:
                isotope = ''
            parts.append(f"{sign}{isotope}i")
        if self.charge != 0 and self.charge != 1:
            charge = abs(self.charge)
            parts.append(f"^{charge}")
        if self.adducts:
            parts.append('[{}]'.format(combine_formula(self.adducts, leading_sign=False)))
        if self.mass_error is not None:
            parts.append("/")
            parts.append(self.mass_error.serialize())
        if self.confidence is not None:
            parts.append(f"*{self.confidence}")
        if self.rest is not None:
            parts.append("/")
            parts.append(self.rest)
        result = ''.join(parts)
        if self.is_auxiliary:
            return f'&{result}'
        return result

    def __str__(self):
        return self.serialize()

    def _molecule_description(self) -> JSONDict:
        return {
            'series_label': self.series_label
        }

    def to_json(self, exclude_missing=False) -> JSONDict:
        #TODO: When neutral losses and adducts are formalized types, convert to string/JSON here
        d = {}
        skips = ('series', 'rest', 'is_auxiliary')
        for key in IonAnnotationBase.__slots__:
            if key in skips:
                continue
            if key == 'mass_error' and self.mass_error is not None:
                d[key] = self.mass_error.to_json()
            else:
                value = getattr(self, key)
                if (value is not None) or not exclude_missing:
                    d[key] = value
        d['molecule_description'] = self._molecule_description()
        if d['analyte_reference'] is None:
            d['analyte_reference'] = '1'
        return d

    def _populate_from_dict(self, data) -> 'IonAnnotationBase':
        #TODO: When neutral losses and adducts are formalized types, parse from string here
        for key, value in data.items():
            if key == 'molecule_description':
                continue
            elif key == 'mass_error' and value is not None:
                self.mass_error = MassError(value['value'], value['unit'])
            else:
                setattr(self, key, value)
        self.rest = None
        self.is_auxiliary = False
        return self

    @classmethod
    def from_json(cls, data) -> 'IonAnnotationBase':
        descr = data["molecule_description"]
        series_label = descr['series_label']
        cls = cls._label_registry[series_label]
        self = cls.__new__(cls)
        self._populate_from_dict(data)
        return self


class PeptideFragmentIonAnnotation(IonAnnotationBase):
    __slots__ = ("position", )

    series_label = 'peptide'

    _molecule_description_fields = {
        "series": "The peptide ion series this ion belongs to",
        "position": "The position from the appropriate terminal along the peptide this ion was fragmented at (starting with 1)"
    }

    position: int

    def __init__(self, series, position, neutral_losses=None, isotope=None, adducts=None, charge=None,
                 analyte_reference=None, mass_error=None, confidence=None, rest=None, is_auxiliary=None):
        super(PeptideFragmentIonAnnotation, self).__init__(
            series, neutral_losses, isotope, adducts, charge, analyte_reference, mass_error, confidence,
            rest, is_auxiliary)
        self.position = position

    def _format_ion(self) -> str:
        return f"{self.series}{self.position}"

    def _molecule_description(self) -> JSONDict:
        d = super()._molecule_description()
        d.update({
            "series": self.series,
            "position": self.position
        })
        return d

    def _populate_from_dict(self, data) -> IonAnnotationBase:
        super()._populate_from_dict(data)
        descr = data['molecule_description']
        self.series = descr['series']
        self.position = descr['position']
        return self


class InternalPeptideFragmentIonAnnotation(IonAnnotationBase):
    __slots__ = ("start_position", "end_position")

    series_label = 'internal'

    _molecule_description_fields = {
        "start_position": ("N-terminal amino acid residue of the fragment in the "
                           "original peptide sequence (beginning with 1, counting from the N-terminus)"),
        "end_position": ("C-terminal amino acid residue of the fragment in the original peptide sequence "
                         "(beginning with 1, counting from the N-terminus)")
    }

    start_position: int
    end_position: int

    def __init__(self, series, start_position, end_position, neutral_losses=None, isotope=None,
                 adducts=None, charge=None, analyte_reference=None, mass_error=None, confidence=None,
                 rest=None, is_auxiliary=None):
        super(InternalPeptideFragmentIonAnnotation, self).__init__(
            series, neutral_losses, isotope, adducts, charge, analyte_reference, mass_error, confidence,
            rest, is_auxiliary)
        self.start_position = start_position
        self.end_position = end_position

    def _format_ion(self) -> str:
        return f"m{self.start_position}:{self.end_position}"

    def _molecule_description(self) -> JSONDict:
        d = super()._molecule_description()
        d['start_position'] = self.start_position
        d['end_position'] = self.end_position
        return d

    def _populate_from_dict(self, data) -> IonAnnotationBase:
        super()._populate_from_dict(data)
        descr = data['molecule_description']
        self.start_position = descr['start_position']
        self.end_position = descr['end_position']
        return self


class PrecursorIonAnnotation(IonAnnotationBase):
    __slots__ = ()

    series_label = "precursor"
    _molecule_description_fields = {}

    def __init__(self, series, neutral_losses=None, isotope=None, adducts=None, charge=None,
                 analyte_reference=None, mass_error=None, confidence=None, rest=None,
                 is_auxiliary=None):
        super(PrecursorIonAnnotation, self).__init__(
            series, neutral_losses, isotope, adducts, charge, analyte_reference, mass_error,
            confidence, rest, is_auxiliary)

    def _format_ion(self):
        return "p"


class ImmoniumIonAnnotation(IonAnnotationBase):
    __slots__ = ("amino_acid", "modification")

    series_label = "immonium"
    _molecule_description_fields = {
        "amino_acid": "The amino acid represented by this immonium ion",
        "modification": "An optional modification that may be attached to this immonium ion"
    }

    amino_acid: str
    modification: str

    def __init__(self, series, amino_acid, modification=None, neutral_losses=None, isotope=None, adducts=None,
                 charge=None, analyte_reference=None, mass_error=None, confidence=None, rest=None, is_auxiliary=None):
        super(ImmoniumIonAnnotation, self).__init__(
            series, neutral_losses, isotope, adducts, charge, analyte_reference, mass_error, confidence,
            rest, is_auxiliary)
        self.amino_acid = amino_acid
        self.modification = modification

    def _format_ion(self):
        if self.modification is not None:
            modification = f"[{self.modification}]"
        else:
            modification = ''
        return f"I{self.amino_acid}{modification}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d['amino_acid'] = self.amino_acid
        if self.modification:
            d['modification'] = self.modification
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data['molecule_description']
        self.amino_acid = descr['amino_acid']
        self.modification = descr.get('modification')
        return self


class ReporterIonAnnotation(IonAnnotationBase):
    __slots__ = ("reporter_label", )

    series_label = "reporter"
    _molecule_description_fields = {
        "reporter_label": "The labeling reagent's name or channel information"
    }

    reporter_label: str

    def __init__(self, series, reporter_label, neutral_losses=None, isotope=None, adducts=None, charge=None,
                 analyte_reference=None, mass_error=None, confidence=None, rest=None, is_auxiliary=None):
        super(ReporterIonAnnotation, self).__init__(
            series, neutral_losses, isotope, adducts, charge, analyte_reference, mass_error, confidence,
            rest, is_auxiliary)
        self.reporter_label = reporter_label

    def _format_ion(self):
        return f"r[{self.reporter_label}]"

    def _molecule_description(self):
        d = super()._molecule_description()
        d['reporter_label'] = self.reporter_label
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data['molecule_description']
        self.reporter_label = descr['reporter_label']
        return self


class ExternalIonAnnotation(IonAnnotationBase):
    __slots__ = ('label', )

    series_label = "external"

    _molecule_description_fields = {
        "label": "The name of the external ion being marked"
    }

    label: str

    def __init__(self, series, label, neutral_losses=None, isotope=None, adducts=None, charge=None,
                 analyte_reference=None, mass_error=None, confidence=None, rest=None, is_auxiliary=None):
        super(ExternalIonAnnotation, self).__init__(
            series, neutral_losses, isotope, adducts, charge, analyte_reference, mass_error, confidence,
            rest, is_auxiliary)
        self.label = label

    def _format_ion(self):
        return f"_{{{self.label}}}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d['label'] = self.label
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data['molecule_description']
        self.label = descr['label']
        return self


class FormulaAnnotation(IonAnnotationBase):
    __slots__ = ("formula", )

    series_label = "formula"
    _molecule_description_fields = {
        "formula": "The elemental formula of the ion being marked"
    }

    formula: str

    def __init__(self, series, formula, neutral_losses=None, isotope=None, adducts=None, charge=None,
                 analyte_reference=None, mass_error=None, confidence=None, rest=None, is_auxiliary=None):
        super(FormulaAnnotation, self).__init__(
            series, neutral_losses, isotope, adducts, charge, analyte_reference, mass_error, confidence,
            rest, is_auxiliary)
        self.formula = formula

    def _format_ion(self):
        return f"f{{{self.formula}}}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d['formula'] = self.formula
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data['molecule_description']
        self.formula = descr['formula']
        return self


class SMILESAnnotation(IonAnnotationBase):
    __slots__ = ("smiles", )

    series_label = "smiles"
    _molecule_description_fields = {
        "smiles": "The SMILES definition of the ion being marked"
    }

    smiles: str

    def __init__(self, series, smiles, neutral_losses=None, isotope=None, adducts=None, charge=None,
                 analyte_reference=None, mass_error=None, confidence=None, rest=None, is_auxiliary=None):
        super().__init__(
            series, neutral_losses, isotope, adducts, charge, analyte_reference, mass_error, confidence,
            rest, is_auxiliary)
        self.smiles = smiles

    def _format_ion(self):
        return f"s{{{self.smiles}}}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d['smiles'] = self.smiles
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data['molecule_description']
        self.smiles = descr['smiles']
        return self


class Unannotated(IonAnnotationBase):
    series_label = "unannotated"
    _molecule_description_fields = {
        "unannotated_label": "A user-specified numeral label for an unannotated peak"
    }

    __slots__ = ('unannotated_label', )

    unannotated_label: str

    def __init__(self, series, unannotated_label, neutral_losses=None, isotope=None, adducts=None, charge=None,
                 analyte_reference=None, mass_error=None, confidence=None, rest=None, is_auxiliary=None):
        self.unannotated_label = unannotated_label
        super().__init__(
            series, neutral_losses, isotope, adducts, charge, analyte_reference,
            mass_error, confidence, rest, is_auxiliary)

    def _format_ion(self):
        if not self.unannotated_label:
            return "?"
        return f"?{self.unannotated_label}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d['unannotated_label'] = self.unannotated_label
        return d

    def serialize(self) -> str:
        # mass_error is a required field in the data model, but it is meaningless for unannotated peaks,
        # so we mask it here
        mass_error = self.mass_error

        mask = mass_error.mass_error == 0
        if mask:
            self.mass_error = None

        val = super().serialize()

        if mask:
            self.mass_error = mass_error
        return val

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data['molecule_description']
        self.unannotated_label = descr.get('unannotated_label')
        return self


class InvalidAnnotation(IonAnnotationBase):
    series_label = "!invalid!"

    content: str
    error: str

    def __init__(self, content: str, error: str):
        super().__init__(self.series_label)
        self.content = content
        self.error = error

    def serialize(self):
        return self.content

    def _molecule_description(self) -> JSONDict:
        return {
            "series_label": self.series_label,
            "content": self.content,
            "error": self.error
        }

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data['molecule_description']
        self.content = descr['content']
        self.error = descr['error']
        return self


def int_or_sign(string: str) -> int:
    if string == "+":
        return 1
    elif string == '-':
        return -1
    else:
        return int(string)


class AnnotationStringParser(object):
    pattern: Pattern

    def __init__(self, pattern):
        self.pattern = pattern

    def __call__(self, annotation_string: str, *, wrap_errors=True, **kwargs) -> List[IonAnnotationBase]:
        try:
            return self.parse_annotation(annotation_string, **kwargs)
        except ValueError as err:
            return [InvalidAnnotation(annotation_string, str(err))]

    def _parse_string(self, annotation_string: str, **kwargs) -> Tuple[re.Match, Dict[str, str]]:
        match = self.pattern.search(annotation_string)
        if match is None:
            raise ValueError(f"Invalid annotation string {annotation_string!r}")
        data = match.groupdict()
        return match, data

    def _coerce_isotope(self, data: Dict[str, str]) -> int:
        isotope = int_or_sign(data.get('isotope', 0) or 0)
        return isotope

    def _coerce_charge(self, data: Dict[str, str]) -> int:
        charge = (data.get("charge", 1))
        if charge is None:
            charge = 1
        elif charge == 0:
            raise ValueError(
                f"The charge of an annotation cannot be zero (parsed {data['charge']})")
        else:
            charge = int(charge)
        return charge

    def _coerce_adducts(self, data: Dict[str, str]) -> List[str]:
        adducts = tokenize_signed_symbol_list(data.get("adducts"))
        return adducts

    def _coerce_analyte_reference(self, data: Dict[str, str]) -> str:
        return data.get("analyte_reference", '1')

    def _coerce_neutral_losses(self, data: Dict[str, str]) -> List:
        return tokenize_signed_symbol_list(data.get("neutral_losses"))

    def _coerce_mass_error(self, data: Dict[str, str]) -> MassError:
        mass_error = data.get("mass_error")
        if mass_error is not None:
            mass_error = MassError(float(mass_error), data.get("mass_error_unit"))
        return mass_error

    def _coerce_confidence(self, data: Dict[str, str]) -> float:
        confidence = data.get('confidence')
        if confidence is not None:
            confidence = float(confidence)
            if confidence > 1.0:
                raise ValueError(
                    "A single peak interpretation's confidence cannot be greater than 1.")
        return confidence

    def parse_annotation(self, annotation_string: str, **kwargs) -> List[IonAnnotationBase]:
        if not annotation_string:
            return []

        match, data = self._parse_string(annotation_string)

        is_auxiliary = bool(data.get('is_auxiliary'))
        adducts = self._coerce_adducts(data)
        charge = self._coerce_charge(data)
        isotope = self._coerce_isotope(data)

        # FIXME: ensure that neutral loss is not a plain mass, and tokenize separate blocks
        neutral_losses = self._coerce_neutral_losses(data)

        analyte_reference =self._coerce_analyte_reference(data)
        mass_error = self._coerce_mass_error(data)
        confidence = self._coerce_confidence(data)

        annotation = self._dispatch(
            annotation_string,
            data,
            adducts,
            charge,
            isotope,
            neutral_losses,
            analyte_reference,
            mass_error,
            confidence,
            **kwargs
        )
        if is_auxiliary:
            annotation.is_auxiliary = True

        rest = annotation_string[match.end():]
        if rest == "":
            return [annotation]
        else:
            if rest[0] != ",":
                raise ValueError(f"Malformed trailing string {rest}, expected ',' for {annotation_string}")
            else:
                rest = rest[1:]
            result = [annotation]
            result.extend(self.parse_annotation(rest, **kwargs))
            total_confidence = 0.0
            for annot in result:
                if annot.confidence is not None:
                    total_confidence += annot.confidence
            if total_confidence < 0 or total_confidence > (1 + 1e-3):
                raise ValueError(
                    f"The sum of all interpretations of a single peak's confidence cannot be greater than 1"
                    f" ({total_confidence}). {annotation_string}")
            return result

    def _dispatch(self, annotation_string, data, adducts, charge, isotope, neutral_losses, analyte_reference,
                  mass_error, confidence, **kwargs):
        if data.get('series'):
            return self._dispatch_peptide_fragment(
                data,
                neutral_losses=neutral_losses, isotope=isotope, adducts=adducts, charge=charge,
                analyte_reference=analyte_reference, mass_error=mass_error, confidence=confidence, **kwargs)
        elif data.get('series_internal'):
            return self._dispatch_internal_peptide_fragment(
                data,
                neutral_losses=neutral_losses, isotope=isotope, adducts=adducts, charge=charge,
                analyte_reference=analyte_reference, mass_error=mass_error, confidence=confidence, **kwargs)
        elif data.get('precursor'):
            return self._dispatch_precursor(
                data,
                neutral_losses=neutral_losses, isotope=isotope, adducts=adducts, charge=charge,
                analyte_reference=analyte_reference, mass_error=mass_error, confidence=confidence, **kwargs)
        elif data.get('immonium'):
            return self._dispatch_immonium(
                data,
                neutral_losses=neutral_losses, isotope=isotope, adducts=adducts, charge=charge,
                analyte_reference=analyte_reference, mass_error=mass_error, confidence=confidence, **kwargs)
        elif data.get('reporter'):
            return self._dispatch_reporter(
                data,
                neutral_losses=neutral_losses, isotope=isotope, adducts=adducts, charge=charge,
                analyte_reference=analyte_reference, mass_error=mass_error, confidence=confidence, **kwargs)
        elif data.get('external_ion'):
            return self._dispatch_external(
                data,
                neutral_losses=neutral_losses, isotope=isotope, adducts=adducts, charge=charge,
                analyte_reference=analyte_reference, mass_error=mass_error, confidence=confidence, **kwargs)
        elif data.get('formula'):
            return self._dispatch_formula(
                data,
                neutral_losses=neutral_losses, isotope=isotope, adducts=adducts, charge=charge,
                analyte_reference=analyte_reference, mass_error=mass_error, confidence=confidence, **kwargs
            )
        elif data.get("smiles"):
            return self._dispatch_smiles(
                data,
                neutral_losses=neutral_losses, isotope=isotope, adducts=adducts, charge=charge,
                analyte_reference=analyte_reference, mass_error=mass_error, confidence=confidence, **kwargs
            )
        elif data.get('unannotated'):
            return self._dispatch_unannotated(
                data,
                neutral_losses=neutral_losses, isotope=isotope, adducts=adducts, charge=charge,
                analyte_reference=analyte_reference, mass_error=mass_error, confidence=confidence, **kwargs
            )
        else:
            raise ValueError(f"Could not infer annotation type from {annotation_string}/{data}")

    def _dispatch_peptide_fragment(self, data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, confidence, **kwargs):
        return PeptideFragmentIonAnnotation(
            data['series'], int(data['ordinal']),
            neutral_losses, isotope, adducts, charge, analyte_reference,
            mass_error, confidence)

    def _dispatch_unannotated(self, data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, confidence, **kwargs):
        if mass_error is None:
            mass_error = MassError(0)
        return Unannotated(
            None, data.get('unannotated_label'), neutral_losses, isotope, adducts, charge, analyte_reference,
            mass_error, confidence)

    def _dispatch_internal_peptide_fragment(self, data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, confidence, **kwargs):
        return InternalPeptideFragmentIonAnnotation(
            "internal", int(data['internal_start']), int(data['internal_end']),
            neutral_losses, isotope, adducts, charge, analyte_reference,
            mass_error, confidence)

    def _dispatch_precursor(self, data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, confidence, **kwargs):
        return PrecursorIonAnnotation(
            "precursor",
            neutral_losses, isotope, adducts, charge, analyte_reference,
            mass_error, confidence)

    def _dispatch_immonium(self, data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, confidence, **kwargs):
        return ImmoniumIonAnnotation(
            "immonium", data['immonium'], data['immonium_modification'],
            neutral_losses, isotope, adducts, charge, analyte_reference,
            mass_error, confidence)

    def _dispatch_reporter(self, data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, confidence, **kwargs):
        return ReporterIonAnnotation(
            "reporter", (data["reporter_label"]),
            neutral_losses, isotope, adducts, charge, analyte_reference,
            mass_error, confidence)

    def _dispatch_external(self, data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, confidence, **kwargs):
        return ExternalIonAnnotation(
            "external", data['external_ion'],
            neutral_losses, isotope, adducts, charge, analyte_reference,
            mass_error, confidence)

    def _dispatch_formula(self, data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, confidence, **kwargs):
            return FormulaAnnotation(
                "formula", data['formula'],
                neutral_losses, isotope, adducts, charge, analyte_reference,
                mass_error, confidence)

    def _dispatch_smiles(self, data, adducts, charge, isotope, neutral_losses, analyte_reference, mass_error, confidence, **kwargs):
            return SMILESAnnotation(
                "smiles", data['smiles'],
                neutral_losses, isotope, adducts, charge, analyte_reference,
                mass_error, confidence)

parse_annotation = AnnotationStringParser(annotation_pattern)


if __name__ == "__main__":
    print(_sre_to_ecma(annotation_pattern.pattern))
