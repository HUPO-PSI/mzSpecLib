import io
import re

from typing import Dict, Tuple

from mzlib.attributes import AttributeManager

from mzlib.annotation import AnnotationStringParser

from mzlib.backends.base import DEFAULT_VERSION, FORMAT_VERSION_TERM, LIBRARY_NAME_TERM
from mzlib.backends.utils import open_stream

from .msp import MSPSpectralLibrary as _MSPSpectralLibrary


annotation_pattern = re.compile(r"""^
(?:(?:(?P<series>[abyxcz]\.?)(?P<ordinal>\d+))|
   (:?Int/(?P<series_internal>[ARNDCEQGHKMFPSTWYVILJarndceqghkmfpstwyvilj]+))|
   (?P<precursor>p)|
   (:?I(?P<immonium>[ARNDCEQGHKMFPSTWYVIL])(?:(?P<immonium_modification>CAM)|[A-Z])?)|
   (?:_(?P<external_ion>[^\s,/]+))
)
(?P<neutral_losses>[+-]\d*(?:\.\d+)?)
(?:\^(?P<charge>[+-]?\d+))?
(?:\[M(?P<adducts>(:?[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?
(?:(?P<isotope>\d*)i)?
(?:/(?P<mass_error>[+-]?\d+(?:\.\d+))(?P<mass_error_unit>ppm)?)?
""", re.X)


class SPTXTAnnotationParser(AnnotationStringParser):
    def __init__(self, pattern):
        super().__init__(pattern)

    def _parse_string(self, annotation_string: str, **kwargs) -> Tuple[re.Match, Dict[str, str]]:
        match, data = super()._parse_string(annotation_string, **kwargs)
        if data['isotope'] is not None and not data['isotope']:
            data['isotope'] = '1'
        return match, data


parse_annotation = SPTXTAnnotationParser(annotation_pattern)


sptxt_spectrum_attribute_map = {
    "TotalIonCurrent": "MS:1000285|total ion current",

}


class SPTXTSpectralLibrary(_MSPSpectralLibrary):
    file_format = "sptxt"
    format_name = "sptxt"

    @classmethod
    def guess_from_header(cls, filename: str) -> bool:
        with open_stream(filename, 'r') as stream:
            first_line = stream.readline()
            i = 1
            while first_line.startswith("###") and i < 1000:
                first_line = stream.readline()
                i += 1
            if re.match("Name: ", first_line):
                return True
            if first_line.startswith("###"):
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

    def _parse_annotation(self, annotation: str, wrap_errors: bool = True, **kwargs):
        return parse_annotation(annotation_string=annotation, wrap_errors=wrap_errors, **kwargs)
