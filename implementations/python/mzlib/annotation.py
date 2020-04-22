import re

annotation_pattern = re.compile(r"""
^(?:(?:(?P<series>[bycz]\.?)(?P<ordinal>\d+))|
   (?P<series_internal>[m](?P<internal_start>\d+):(?P<internal_end>\d+))|
   (?P<precursor>p)|
   (?P<immonium>[ARNDCEQGHKMFPSTWYV]|IL|LI|I|L)|
   (?P<reporter>r(?P<reporter_mass>\d+(?:\.\d+)))|
   (?:_(?P<external_ion>[^\s,/]+))
)
(?P<neutral_loss>[+-](?:[A-Z][A-Za-z0-9]+))?
(?:(?P<isotope>[+-]\d*)i)?
(?:\^(?P<charge>[+-]?\d+))?
(?:@(?P<analyte_reference>[^/\s]+))?
(?:/(?P<mass_error>-?\d+(?:\.\d+))(?P<mass_error_unit>ppm)?)?
""", re.X)

# At the time of first writing, this pattern could be translated into the equivalent
# ECMAScript compliant regex:
# (?:(?:(?<series>[bycz])(?<ordinal>\d+))|(?<series_internal>[m](?<internal_start>\d+):(?<internal_end>\d+))|(?<precursor>p)|(?<immonium>[ARNDCEQGHKMFPSTWYV]|IL|LI|I|L)|(?<reporter>r(?<reporter_mass>\d+(?:\.\d+)))|(?:_(?<external_ion>[^\s,/]+)))(?<neutral_loss>[+-](?:[A-Z][A-Za-z0-9]+))?(?:(?<isotope>[+-]\d*i))?(?:\^(?<charge>[+-]?\d+))?(?:@(?<analyte_reference>[^/\s]+))?(?:/(?<mass_error>-?\d+(?:\.\d+))(?<mass_error_unit>ppm)?)?
# Line breaks not introduced to preserve syntactic correctness.


class MassError(object):
    _DEFAULT_UNIT = "Da"

    def __init__(self, mass_error, unit=None):
        if unit is None:
            unit = self._DEFAULT_UNIT
        self.mass_error = float(mass_error)
        self.unit = unit

    def serialize(self):
        unit = self.unit
        if unit == self._DEFAULT_UNIT:
            unit = ''
        return f"{self.mass_error}{unit}"

    def __str__(self):
        return self.serialize()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.mass_error}, {self.unit})"


class IonAnnotationBase(object):
    __slots__ = ("series", "neutral_loss", "isotope", "adduct", "charge", "analyte_reference",
                 "mass_error", "rest")

    def __init__(self, series, neutral_loss=None, isotope=None, adduct=None, charge=None,
                 analyte_reference=None,
                 mass_error=None, rest=None):
       if isotope is None:
            isotope = 0
       if charge is None:
            charge = 1
       self.series = series
       self.neutral_loss = neutral_loss
       self.isotope = isotope
       self.adduct = adduct
       self.charge = charge
       self.analyte_reference = analyte_reference
       self.mass_error = mass_error
       self.rest = rest

    def __hash__(self):
        return hash(self.serialize())

    def __eq__(self, other):
        return self.serialize() == str(other)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return self.serialize()

    def _format_ion(self):
        raise NotImplementedError()

    def serialize(self):
        parts = [self._format_ion(), ]
        if self.neutral_loss is not None:
            parts.append(str(self.neutral_loss))
        if self.isotope != 0:
            sign = "+" if self.isotope > 0 else "-"
            isotope = abs(self.isotope)
            if isotope == 1:
                isotope = ''
            parts.append(f"{sign}{isotope}i")
        if self.adduct is not None:
            parts.append(self.adduct)
        if self.charge != 0 and self.charge != 1:
            charge = abs(self.charge)
            parts.append(f"^{charge}")
        if self.analyte_reference is not None:
            parts.append(f"@{self.analyte_reference}")
        if self.mass_error is not None:
            parts.append("/")
            parts.append(self.mass_error.serialize())
        if self.rest is not None:
            parts.append("/")
            parts.append(self.rest)
        return ''.join(parts)

    def __str__(self):
        return self.serialize()


class PeptideFragmentIonAnnotaiton(IonAnnotationBase):
    def __init__(self, series, position, neutral_loss=None, isotope=None, adduct=None, charge=None,
                 analyte_reference=None, mass_error=None, rest=None):
        super(PeptideFragmentIonAnnotaiton, self).__init__(
            series, neutral_loss, isotope, adduct, charge, analyte_reference, mass_error, rest)
        self.position = position

    def _format_ion(self):
        return f"{self.series}{self.position}"


class InternalPeptideFragmentIonAnnotation(IonAnnotationBase):
    def __init__(self, series, start_position, end_position, neutral_loss=None, isotope=None,
                 adduct=None, charge=None, analyte_reference=None, mass_error=None, rest=None):
        super(InternalPeptideFragmentIonAnnotation, self).__init__(
            series, neutral_loss, isotope, adduct, charge, analyte_reference, mass_error, rest)
        self.start_position = start_position
        self.end_position = end_position

    def _format_ion(self):
        return f"{self.series}{self.start_position}:{self.end_position}"


class PrecursorIonAnnotation(IonAnnotationBase):
    def __init__(self, series, neutral_loss=None, isotope=None, adduct=None, charge=None,
                 analyte_reference=None, mass_error=None, rest=None):
        super(PrecursorIonAnnotation, self).__init__(
            series, neutral_loss, isotope, adduct, charge, analyte_reference, mass_error, rest)

    def _format_ion(self):
        return "p"


class ImmoniumIonAnnotation(IonAnnotationBase):
    def __init__(self, series, amino_acids, neutral_loss=None, isotope=None, adduct=None, charge=None,
                 analyte_reference=None, mass_error=None, rest=None):
        super(ImmoniumIonAnnotation, self).__init__(
            series, neutral_loss, isotope, adduct, charge, analyte_reference, mass_error, rest)
        self.amino_acids = amino_acids

    def _format_ion(self):
        return f"{self.amino_acids}"


class ReporterIonAnnotation(IonAnnotationBase):
    def __init__(self, series, reporter_mass, neutral_loss=None, isotope=None, adduct=None, charge=None,
                 analyte_reference=None, mass_error=None, rest=None):
        super(ReporterIonAnnotation, self).__init__(
            series, neutral_loss, isotope, adduct, charge, analyte_reference, mass_error, rest)
        self.reporter_mass = reporter_mass

    def _format_ion(self):
        return f"r{self.reporter_mass}"


class ExternalIonAnnotation(IonAnnotationBase):
    def __init__(self, series, label, neutral_loss=None, isotope=None, adduct=None, charge=None,
                 analyte_reference=None, mass_error=None, rest=None):
        super(ExternalIonAnnotation, self).__init__(
            series, neutral_loss, isotope, adduct, charge, analyte_reference, mass_error, rest)
        self.label = label

    def _format_ion(self):
        return f"_{self.label}"


def int_or_sign(string):
    if string == "+":
        return 1
    elif string == '-':
        return -1
    else:
        return int(string)


def parse_annotation(annotation_string):
    if annotation_string == "?" or not annotation_string:
        return None
    match = annotation_pattern.search(annotation_string)
    if match is None:
        raise ValueError(f"Invalid annotation string {annotation_string!r}")
    data = match.groupdict()

    adduct = None
    charge = (data.get("charge", 1))
    if charge is None:
        charge = 1
    elif charge == 0:
        raise ValueError(f"The charge of an annotation cannot be zero. {annotation_string}")
    else:
        charge = int(charge)
    isotope = int_or_sign(data.get('isotope', 0) or 0)
    neutral_loss = data.get("neutral_loss")
    analyte_reference = data.get("analyte_reference")

    mass_error = data.get("mass_error")
    if mass_error is not None:
        mass_error = MassError(float(mass_error), data.get("mass_error_unit"))
    if data['series']:
        return PeptideFragmentIonAnnotaiton(
            data['series'], int(data['ordinal']),
            neutral_loss, isotope, adduct, charge, analyte_reference,
            mass_error)
    elif data['series_internal']:
        return InternalPeptideFragmentIonAnnotation(
            "internal", int(data['internal_start']), int(data['internal_end']),
            neutral_loss, isotope, adduct, charge, analyte_reference,
            mass_error)
    elif data['precursor']:
        return PrecursorIonAnnotation(
            "precursor",
            neutral_loss, isotope, adduct, charge, analyte_reference,
            mass_error)
    elif data['immonium']:
        return ImmoniumIonAnnotation(
            "immonium", data['immonium'],
            neutral_loss, isotope, adduct, charge, analyte_reference,
            mass_error)
    elif data['reporter']:
        return ReporterIonAnnotation(
            "reporter", float(data["reporter_mass"]),
            neutral_loss, isotope, adduct, charge, analyte_reference,
            mass_error)
    elif data['external_ion']:
        return ExternalIonAnnotation(
            "external", data['external_ion'],
            neutral_loss, isotope, adduct, charge, analyte_reference,
            mass_error)
    else:
        raise ValueError(f"Could not infer annotation type from {annotation_string}")
