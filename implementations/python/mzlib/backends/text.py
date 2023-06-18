import re
import os
import io
import logging
import warnings
import enum

from collections import deque
from typing import ClassVar, List, Optional, Tuple, Union, Iterable

from mzlib.annotation import parse_annotation
from mzlib.spectrum import Spectrum
from mzlib.cluster import SpectrumCluster
from mzlib.attributes import AttributeManager, Attributed, AttributeSet
from mzlib.analyte import Analyte, Interpretation, InterpretationMember

from .base import (
    SpectralLibraryBackendBase,
    _PlainTextSpectralLibraryBackendBase,
    SpectralLibraryWriterBase,
    FORMAT_VERSION_TERM,
    AttributeSetTypes)
from .utils import try_cast, open_stream

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


term_pattern = re.compile(
    r"^(?P<term>(?P<term_accession>\S+:(?:\d|X)+)\|(?P<term_name>[^=]+))")
key_value_term_pattern = re.compile(
    r"^(?P<term>(?P<term_accession>[A-Za-z0-9:.]+:(?:\d|X)+)\|(?P<term_name>[^=]+))=(?P<value>.+)")
grouped_key_value_term_pattern = re.compile(
    r"^\[(?P<group_id>\d+)\](?P<term>(?P<term_accession>\S+:(?:\d|X)+)\|(?P<term_name>[^=]+))=(?P<value>.+)")
float_number = re.compile(
    r"^\d+(.\d+)?")


class _SpectrumParserStateEnum(enum.Enum):
    unknown = 0
    header = 1
    analyte = 2
    interpretation = 3
    interpretation_member = 4
    peaks = 5
    done = 6
    cluster = 7


class _LibraryParserStateEnum(enum.Enum):
    unknown = 0
    header = 1
    attribute_sets = 2
    content = 3


ATTRIBUTE_SET_NAME = "MS:1003212|library attribute set name"
PEAK_ATTRIBUTE = "MS:1003254|peak attribute"

START_OF_SPECTRUM_MARKER = re.compile(r"^<(?:Spectrum)(?:=(.+))?>")
START_OF_INTERPRETATION_MARKER = re.compile(r"^<Interpretation(?:=(.+))>")
START_OF_ANALYTE_MARKER = re.compile(r"^<Analyte(?:=(.+))>")
START_OF_PEAKS_MARKER = re.compile(r"^<Peaks>")
START_OF_LIBRARY_MARKER = re.compile(r"^<mzSpecLib\s+(.+)>")
SPECTRUM_NAME_PRESENT = re.compile(r'MS:1003061\|(?:library )?spectrum name=')
START_OF_INTERPRETATION_MEMBER_MARKER = re.compile(r"<InterpretationMember(?:=(.+))>")
START_OF_ATTRIBUTE_SET = re.compile(
    r"<AttributeSet (Spectrum|Analyte|Interpretation|Cluster)=(.+)>")
START_OF_CLUSTER = re.compile(r"<Cluster(?:=(.+))>")


attribute_set_types = {
    "spectrum": AttributeSetTypes.spectrum,
    "analyte": AttributeSetTypes.analyte,
    "interpretation": AttributeSetTypes.interpretation,
    "cluster": AttributeSetTypes.cluster,
}


class _EntryParser:
    """
    Moves the complexity and state management involved in parsing
    a full entry out of :class:`TextSpectrumLibrary`, allowing it
    to be factored into a bunch of helper methods around a single
    piece of shared stated too granular for the main parser.
    """

    library: 'TextSpectralLibrary'
    state: _SpectrumParserStateEnum
    spectrum: Optional[Spectrum]
    cluster: Optional[SpectrumCluster]
    analyte: Optional[Analyte]
    interpretation: Optional[Interpretation]
    interpretation_member: Optional[InterpretationMember]

    aggregation_types: List[str]
    peak_list: List[Tuple]

    start_line_number: int
    line_number: int = -1

    def __init__(self, library, start_line_number: int, spectrum_index: Optional[int]) -> None:
        self.library = library
        self.start_line_number = start_line_number
        self.spectrum_index = spectrum_index
        self.state = _SpectrumParserStateEnum.header

        self.aggregation_types = None
        self.peak_list = []

        self.spectrum = None
        self.cluster = None
        self.analyte = None
        self.interpretation = None
        self.interpretation_member = None

    def real_line_number_or_nothing(self):
        if self.start_line_number is None:
            return ''
        message = f" on line {self.line_number + self.start_line_number}"
        if self.spectrum_index is not None:
            message += f" in spectrum {self.spectrum_index}"
        message += f" in state {self.state}"
        return message

    def _parse_header(self, line):
        if START_OF_SPECTRUM_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.header
            self.spectrum = self.library._new_spectrum()
            self.spectrum.index = self.spectrum_index
            match = START_OF_SPECTRUM_MARKER.match(line)
            self.spectrum.key = int(match.group(1)) or self.spectrum.index - 1
            return

        elif START_OF_PEAKS_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.peaks
            return

        elif START_OF_INTERPRETATION_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.interpretation
            match = START_OF_INTERPRETATION_MARKER.match(line)
            if self.interpretation is not None:
                self.spectrum.add_interpretation(self.interpretation)
            self.interpretation = self.library._new_interpretation(match.group(1))
            self.spectrum.add_interpretation(self.interpretation)
            self.analyte = None
            return

        elif START_OF_ANALYTE_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.analyte
            match = START_OF_ANALYTE_MARKER.match(line)
            self.analyte = self.library._new_analyte(match.group(1))
            self.spectrum.add_analyte(self.analyte)
            return

        elif START_OF_CLUSTER.match(line):
            self.state = _SpectrumParserStateEnum.cluster
            self.cluster = self.library._new_cluster()
            match = START_OF_CLUSTER.match(line)
            self.cluster.key = int(match.group(1)) or self.cluster.index - 1
            return

        self.library._parse_attribute_into(
            line, self.spectrum, self.real_line_number_or_nothing, self.state)

    def _parse_interpretation(self, line):
        if START_OF_ANALYTE_MARKER.match(line):
            warnings.warn(
                f"An analyte found after an interpretation was encountered, {self.real_line_number_or_nothing()}")
            self.state = _SpectrumParserStateEnum.analyte
            match = START_OF_ANALYTE_MARKER.match(line)
            if self.analyte is not None:
                self.spectrum.add_analyte(self.analyte)
            self.analyte = self.library._new_analyte(match.group(1))
            self.spectrum.add_analyte(self.analyte)
            return
        elif START_OF_INTERPRETATION_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.interpretation
            match = START_OF_INTERPRETATION_MARKER.match(line)
            if self.interpretation is not None:
                self.spectrum.add_interpretation(self.interpretation)
            self.interpretation = self.library._new_interpretation(match.group(1))
            self.spectrum.add_interpretation(self.interpretation)
            self.analyte = None
            return
        elif START_OF_PEAKS_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.peaks
            return
        elif START_OF_INTERPRETATION_MEMBER_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.interpretation_member
            match = START_OF_INTERPRETATION_MEMBER_MARKER.match(line)

            if self.interpretation_member is not None:
                self.interpretation.add_member_interpretation(self.interpretation_member)

            self.interpretation_member = InterpretationMember(match.group(1))
            self.interpretation.add_member_interpretation(self.interpretation_member)
            return

        self.library._parse_attribute_into(
            line, self.interpretation.attributes, self.real_line_number_or_nothing)
        self.library._analyte_interpretation_link(self.spectrum, self.interpretation)

    def _parse_interpretation_member(self, line):
        if START_OF_PEAKS_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.peaks
            self.interpretation_member = None
            self.interpretation = None
            return
        elif START_OF_INTERPRETATION_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.interpretation
            match = START_OF_INTERPRETATION_MARKER.match(line)
            if self.interpretation is not None:
                self.spectrum.add_interpretation(self.interpretation)
            self.interpretation = self.library._new_interpretation(match.group(1))
            self.spectrum.add_interpretation(self.interpretation)
            self.interpretation_member = None
            return
        elif START_OF_INTERPRETATION_MEMBER_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.interpretation_member
            match = START_OF_INTERPRETATION_MEMBER_MARKER.match(line)
            if self.interpretation_member is not None:
                self.interpretation.add_member_interpretation(self.interpretation_member)
            self.interpretation_member = InterpretationMember(match.group(1))
            self.interpretation.add_member_interpretation(self.interpretation_member)
            return

        self.library._parse_attribute_into(
            line, self.interpretation_member, self.real_line_number_or_nothing)

    def _parse_analyte(self, line):
        if START_OF_PEAKS_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.peaks
            if self.analyte is not None:
                self.spectrum.add_analyte(self.analyte)
                self.analyte = None
            return

        elif START_OF_ANALYTE_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.analyte
            match = START_OF_ANALYTE_MARKER.match(line)
            if self.analyte is not None:
                self.spectrum.add_analyte(self.analyte)
            self.analyte = self.library._new_analyte(match.group(1))
            return

        elif START_OF_INTERPRETATION_MARKER.match(line):
            self.state = _SpectrumParserStateEnum.interpretation
            match = START_OF_INTERPRETATION_MARKER.match(line)
            if self.analyte is not None:
                self.spectrum.add_analyte(self.analyte)
                self.analyte = None

            # Somehow we have an in-progress Interpretation that hasn't been cleared yet.
            # This should probably be an error strictly speaking.
            if self.interpretation is not None:
                warnings.warn(
                    f"Interleaved analytes and interpretations detected at {self.real_line_number_or_nothing()}")
                self.spectrum.add_interpretation(self.interpretation)
            self.interpretation = self.library._new_interpretation(match.group(1))
            self.spectrum.add_interpretation(self.interpretation)
            return

        self.library._parse_attribute_into(line, self.analyte, self.real_line_number_or_nothing)

    def _parse_peaks(self, line):
        # TODO: When we know more about how different aggregations are formatted,
        # look that up here once so we remember it and can use it to process the
        # aggregation columns
        if self.aggregation_types is None:
            self.aggregation_types = self.spectrum.peak_aggregations
        match = float_number.match(line)
        if match is not None:
            tokens = line.split("\t")
            n_tokens = len(tokens)
            if n_tokens == 2:
                mz, intensity = tokens
                annotation = parse_annotation("?")
                self.peak_list.append([float(mz), float(intensity), annotation, []])
            elif n_tokens == 3:
                mz, intensity, annotation = tokens
                if not annotation:
                    annotation = "?"
                annotation = parse_annotation(annotation)
                self.peak_list.append([float(mz), float(intensity), annotation, []])
            elif n_tokens > 3:
                mz, intensity, annotation, *aggregation = tokens
                if not annotation:
                    annotation = "?"
                annotation = parse_annotation(annotation)
                self.peak_list.append(
                    [float(mz), float(intensity), annotation, [try_cast(agg) for agg in aggregation]])
            else:
                raise ValueError(
                    f"Malformed peak line {line} with {n_tokens} entries{self.real_line_number_or_nothing()}")
        else:
            raise ValueError(f"Malformed peak line {line}{self.real_line_number_or_nothing()}")

    def _parse_cluster(self, line):
        if START_OF_SPECTRUM_MARKER.match(line):
            raise ValueError(
                f"Clusters should not include spectrum sections {self.real_line_number_or_nothing()}")

        elif START_OF_PEAKS_MARKER.match(line):
            raise ValueError(
                f"Clusters should not include peaks {self.real_line_number_or_nothing()}")

        elif START_OF_INTERPRETATION_MARKER.match(line):
            raise ValueError(
                f"Clusters should not include interpretation sections {self.real_line_number_or_nothing()}")

        elif START_OF_ANALYTE_MARKER.match(line):
            raise ValueError(
                f"Clusters should not include analyte sections {self.real_line_number_or_nothing()}")

        elif START_OF_INTERPRETATION_MEMBER_MARKER.match(line):
            raise ValueError(
                f"Clusters should not include interpretation member sections {self.real_line_number_or_nothing()}")

        self.library._parse_attribute_into(
            line, self.cluster, self.real_line_number_or_nothing, self.state)

    def parse(self, buffer: Iterable[str]):
        line: str
        for line_number, line in enumerate(buffer):
            self.line_number = line_number
            line = line.strip()
            if not line:
                break
            # Skip comments for now, no round-trip
            if line.startswith("#"):
                continue
            elif self.state == _SpectrumParserStateEnum.header:
                self._parse_header(line)
            elif self.state == _SpectrumParserStateEnum.interpretation:
                self._parse_interpretation(line)
            elif self.state == _SpectrumParserStateEnum.interpretation_member:
                self._parse_interpretation_member(line)
            elif self.state == _SpectrumParserStateEnum.analyte:
                self._parse_analyte(line)
            elif self.state == _SpectrumParserStateEnum.peaks:
                self._parse_peaks(line)
            elif self.state == _SpectrumParserStateEnum.cluster:
                self._parse_cluster(line)
            else:
                raise ValueError(
                    f"Unknown state {self.state}{self.real_line_number_or_nothing()}")
        if self.cluster:
            return self.cluster
        self.spectrum.peak_list = self.peak_list
        # Backfill analytes into interpretations that never explicitly listed them.
        self.library._default_interpretation_to_analytes(self.spectrum)
        return self.spectrum


def _is_header_line(line: str) -> bool:
    if START_OF_SPECTRUM_MARKER.match(line):
        return False
    if START_OF_CLUSTER.match(line):
        return False
    if SPECTRUM_NAME_PRESENT.match(line):
        return False
    return True


class TextSpectralLibrary(_PlainTextSpectralLibraryBackendBase):
    file_format: ClassVar[str] = "mzlb.txt"
    format_name: ClassVar[str] = "text"

    @classmethod
    def guess_from_header(cls, filename: str) -> bool:
        with open_stream(filename, 'r', encoding='utf8') as stream:
            first_line = stream.readline()
            if (START_OF_SPECTRUM_MARKER.match(first_line) or
                START_OF_LIBRARY_MARKER.match(first_line) or
                START_OF_CLUSTER.match(first_line)):
                return True
        return False

    def _parse_header_from_stream(self, stream: io.TextIOBase) -> Tuple[bool, int]:
        nbytes = 0
        first_line = stream.readline()
        nbytes += len(first_line)

        state = _LibraryParserStateEnum.unknown

        current_attribute_set = None
        current_attribute_set_type = None

        if not _is_header_line(first_line):
            return True, 0
        elif START_OF_LIBRARY_MARKER.match(first_line):
            state = _LibraryParserStateEnum.header
            match = START_OF_LIBRARY_MARKER.match(first_line)
            version = match.group(1)
            attributes = AttributeManager()
            attributes.add_attribute(FORMAT_VERSION_TERM, version)
            line = stream.readline()
            while _is_header_line(line):
                nbytes += len(line)
                line = line.strip()
                if not line:
                    line = stream.readline()
                    continue
                match = START_OF_ATTRIBUTE_SET.match(line)
                if match:
                    state = _LibraryParserStateEnum.attribute_sets
                    if current_attribute_set is not None:
                        self._add_attribute_set(
                            current_attribute_set, current_attribute_set_type)

                    current_attribute_set_type = attribute_set_types[
                        match.group(1).lower()]
                    attrib_set_name = match.group(2)
                    current_attribute_set = AttributeSet(attrib_set_name, [])
                else:
                    match = key_value_term_pattern.match(line)

                    # We found a single attribute
                    if match is not None:
                        d = match.groupdict()
                        # If we're in an attribute set, store it in the attribute set
                        if state == _LibraryParserStateEnum.attribute_sets:
                            current_attribute_set.add_attribute(
                                d['term'], try_cast(d['value']))
                        else: # Otherwise store it in the library level attributes
                            attributes.add_attribute(
                                d['term'], try_cast(d['value']))

                        line = stream.readline()
                        nbytes += len(line)
                        continue

                    if line.startswith("["):
                        # We found a grouped attribute
                        match = grouped_key_value_term_pattern.match(line)
                        if match is not None:
                            d = match.groupdict()
                            # If we're in an attribute set, store it in the attribute
                            # set
                            if state == _LibraryParserStateEnum.attribute_sets:
                                current_attribute_set.add_attribute(
                                    d['term'], try_cast(d['value']), d['group_id'])
                                current_attribute_set.group_counter = int(d['group_id'])
                            else:  # Otherwise store it in the library level attributes
                                attributes.add_attribute(
                                    d['term'], try_cast(d['value']), d['group_id'])
                                attributes.group_counter = int(d['group_id'])

                            line = stream.readline()
                            nbytes += len(line)
                            continue
                        else:
                            raise ValueError(
                                f"Malformed grouped attribute {line}")
                    elif "=" in line:
                        name, value = line.split("=", 1)
                        if state == _LibraryParserStateEnum.attribute_sets:
                            current_attribute_set.add_attribute(name, value)
                        else:
                            attributes.add_attribute(name, value)
                    else:
                        raise ValueError(f"Malformed attribute line {line}")
                line = stream.readline()

            if current_attribute_set is not None:
                self._add_attribute_set(
                    current_attribute_set, current_attribute_set_type)
            self.attributes.clear()
            self.attributes._from_iterable(attributes)
            return True, nbytes
        return False, 0

    def read_header(self) -> bool:
        if isinstance(self.filename, io.IOBase):
            return self._parse_header_from_stream(self.filename)[0]
        with open_stream(self.filename, 'rt', encoding='utf8') as stream:
            return self._parse_header_from_stream(stream)[0]

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

        with open_stream(filename, 'rt', encoding='utf8') as infile:
            state = 'header'
            entry_buffer = deque()

            n_spectra = 0
            n_clusters = 0

            start_index = 0
            file_offset = 0

            line_beginning_file_offset = 0
            spectrum_file_offset = 0
            spectrum_name = ''
            current_key = None
            entry_is_cluster = False

            # Required for counting file_offset manually (LF vs CRLF)
            infile.readline()
            file_offset_line_ending = len(infile.newlines) - 1
            infile.seek(0)

            logger.debug(f"Reading {filename} ({file_size} bytes)...")
            while 1:
                line = infile.readline()
                if len(line) == 0:
                    break

                line_beginning_file_offset = file_offset

                #### tell() is twice as slow as counting it myself
                file_offset += len(line) + file_offset_line_ending

                line = line.rstrip()
                if state == 'header':

                    if is_spec := START_OF_SPECTRUM_MARKER.match(line):
                        current_key = int(is_spec.group(1))
                        state = 'body'
                        spectrum_file_offset = line_beginning_file_offset
                        entry_is_cluster = False
                    elif is_clus := START_OF_CLUSTER.match(line):
                        current_key = int(is_clus.group(1))
                        state = 'body'
                        spectrum_file_offset = line_beginning_file_offset
                        entry_is_cluster = True
                    else:
                        continue

                if state == 'body':
                    if len(line) == 0:
                        continue

                    is_spec = START_OF_SPECTRUM_MARKER.match(line)
                    is_clus = START_OF_CLUSTER.match(line)
                    if (is_spec) or (is_clus):
                        if len(entry_buffer) > 0:
                            if not entry_is_cluster:
                                if not spectrum_name:
                                    raise ValueError("No spectrum name")
                                self.index.add(
                                    number=current_key,
                                    offset=spectrum_file_offset,
                                    name=spectrum_name,
                                    analyte=None)
                                n_spectra += 1
                                current_key = int(is_spec.group(1)) if is_spec else int(is_clus.group(1))
                                #### Commit every now and then
                                if n_spectra % 10000 == 0:
                                    self.index.commit()
                                    logger.info(
                                        f"Processed {file_offset} bytes, {n_spectra} spectra read, {n_clusters} read")
                            else:
                                self.index.add_cluster(number=current_key, offset=spectrum_file_offset)
                                if n_clusters % 10000 == 0:
                                    self.index.commit()
                                    logger.info(
                                        f"Processed {file_offset} bytes, {n_spectra} spectra read, {n_clusters} read")
                                n_clusters += 1
                                current_key = int(is_spec.group(1)) if is_spec else int(is_clus.group(1))

                        entry_buffer.clear()
                        entry_is_cluster = bool(is_clus)
                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = ''
                    if re.match(r'MS:1003061\|(?:library )?spectrum name', line):
                        spectrum_name = re.match(r'MS:1003061\|(?:library )?spectrum name=(.+)', line).group(1)

                    entry_buffer.append(line)


            if spectrum_name:
                self.index.add(
                    number=current_key,
                    offset=spectrum_file_offset,
                    name=spectrum_name,
                    analyte=None)
                self.index.commit()
                n_spectra += 1
                logger.info(
                    f"Processed {file_offset} bytes, {n_spectra} spectra read, {n_clusters} read")
            elif entry_is_cluster:
                self.index.add_cluster(
                    number=current_key,
                    offset=spectrum_file_offset,
                )
                self.index.commit()
                n_clusters += 1
                logger.info(
                    f"Processed {file_offset} bytes, {n_spectra} spectra read, {n_clusters} read")

            #### Flush the index
            self.index.commit()

        return n_spectra

    def _buffer_from_stream(self, infile: io.IOBase) -> List:
        state = 'body'
        spectrum_buffer = []

        for line in infile:
            line = line.rstrip()
            if state == 'body':
                if len(line) == 0:
                    continue
                if START_OF_SPECTRUM_MARKER.match(line) or START_OF_CLUSTER.match(line):
                    if len(spectrum_buffer) > 0:
                        return spectrum_buffer
                spectrum_buffer.append(line)
        return spectrum_buffer

    def _prepare_attribute_dict(self, match):
        key = match['term_accession']
        value = match['value']
        try:
            term = self.find_term_for(key)
            match['value'] = term.value_type(value)
        except KeyError:
            match['value'] = try_cast(value)

    def _parse_attribute_into(self, line: str, store: Attributed,
                              line_number_message=lambda:'',
                              state: _SpectrumParserStateEnum=None) -> bool:
        match = key_value_term_pattern.match(line)
        if match is not None:
            d = match.groupdict()
            self._prepare_attribute_dict(d)
            if d['term'] == ATTRIBUTE_SET_NAME:
                if _SpectrumParserStateEnum.header == state:
                    attr_set = self.spectrum_attribute_sets[d['value']]
                elif _SpectrumParserStateEnum.analyte == state:
                    attr_set = self.analyte_attribute_sets[d['value']]
                elif _SpectrumParserStateEnum.interpretation == state:
                    attr_set = self.interpretation_attribute_sets[d['value']]
                elif _SpectrumParserStateEnum.cluster == state:
                    attr_set = self.cluster_attribute_sets[d['value']]
                else:
                    raise ValueError(f"Cannot define attribute sets for {state}")
                attr_set.apply(store)
            else:
                store.add_attribute(d['term'], try_cast(d['value']))
            return True
        if line.startswith("["):
            match = grouped_key_value_term_pattern.match(line)
            if match is not None:
                d = match.groupdict()
                self._prepare_attribute_dict(d)
                store.add_attribute(
                    d['term'], try_cast(d['value']), d['group_id'])
                store.group_counter = int(d['group_id'])
                return True
            else:
                raise ValueError(
                    f"Malformed grouped attribute {line}{line_number_message()}")
        elif "=" in line:
            name, value = line.split("=", 1)
            store.add_attribute(name, try_cast(value))
            return True
        else:
            raise ValueError(f"Malformed attribute line {line}{line_number_message()}")

    def _parse(self, buffer: Iterable[str], spectrum_index: int = None,
               start_line_number: int=None) -> Union[Spectrum, SpectrumCluster]:
        parser = _EntryParser(self, start_line_number, spectrum_index)
        return parser.parse(buffer)

    def get_spectrum(self, spectrum_number: int=None,
                     spectrum_name: str=None) -> Spectrum:
        # keep the two branches separate for the possibility that this is not
        # possible with all index schemes.
        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError(
                    "Provide only one of spectrum_number or spectrum_name")
            index_record = self.index.record_for(spectrum_number)
            offset = index_record.offset
        elif spectrum_name is not None:
            index_record = self.index.record_for(spectrum_name)
            offset = index_record.offset
            spectrum_number = index_record.number

        buffer = self._get_lines_for(offset)
        spectrum = self._parse(buffer, index_record.index)
        return spectrum

    def get_cluster(self, cluster_number: int) -> SpectrumCluster:
        offset = self.index.offset_for_cluster(cluster_number)
        buffer = self._get_lines_for(offset)
        cluster = self._parse(buffer, cluster_number)
        return cluster


class TextSpectralLibraryWriter(SpectralLibraryWriterBase):
    file_format = "mzlb.txt"
    format_name = "text"
    default_version = '1.0'

    def __init__(self, filename, version=None, compact_interpretations: bool=True, **kwargs):
        super(TextSpectralLibraryWriter, self).__init__(filename)
        self.version = version
        self._coerce_handle(self.filename)
        self.compact_interpretations = compact_interpretations

    def _write_attributes(self, attributes: Attributed):
        for attribute in attributes:
            value = attribute.value
            if ":" in attribute.key:
                try:
                    term = self.find_term_for(attribute.key.split("|")[0])
                    value = term.value_type.format(value)
                except (KeyError, ValueError):
                    pass
            if attribute.group_id is None:
                self.handle.write(f"{attribute.key}={value}\n")
            else:
                self.handle.write(
                    f"[{attribute.group_id}]{attribute.key}={value}\n")

    def write_header(self, library: SpectralLibraryBackendBase):
        if self.version is None:
            version = library.attributes.get_by_name("format version")
            if version is None:
                version = self.default_version
        else:
            version = self.version
        self.handle.write("<mzSpecLib %s>\n" % (version, ))
        self._write_attributes(
            self._filter_attributes(library.attributes, lambda x: x.key != FORMAT_VERSION_TERM)
        )
        for attr_set in library.spectrum_attribute_sets.values():
            self.write_attribute_set(attr_set, AttributeSetTypes.spectrum)

        for attr_set in library.analyte_attribute_sets.values():
            self.write_attribute_set(attr_set, AttributeSetTypes.analyte)

        for attr_set in library.interpretation_attribute_sets.values():
            self.write_attribute_set(attr_set, AttributeSetTypes.interpretation)

    def write_attribute_set(self, attribute_set: AttributeSet,
                            attribute_set_type: AttributeSetTypes):
        if attribute_set_type == AttributeSetTypes.spectrum:
            set_type = "Spectrum"
        elif attribute_set_type == AttributeSetTypes.analyte:
            set_type = "Analyte"
        elif attribute_set_type == AttributeSetTypes.interpretation:
            set_type = "Interpretation"
        elif attribute_set_type == AttributeSetTypes.cluster:
            set_type = "Cluster"

        header = f"<AttributeSet {set_type}={attribute_set.name}>\n"
        self.handle.write(header)
        self._write_attributes(attribute_set.attributes)
        if attribute_set.attributes:
            self.handle.write('\n')

    def write_spectrum(self, spectrum: Spectrum):
        self.handle.write(f"<Spectrum={spectrum.key}>\n")
        attribs_of = list(self._filter_attributes(
            spectrum,
            self._not_entry_key_or_index)
        )
        self._write_attributes(attribs_of)
        for analyte in spectrum.analytes.values():
            self.handle.write(f"<Analyte={analyte.id}>\n")
            self._write_attributes(analyte.attributes)
        _n_interps = len(spectrum.interpretations)
        for interpretation in spectrum.interpretations.values():
            interpretation: Interpretation

            if len(spectrum.analytes) == 1:
                attribs_of = list(self._filter_attributes(
                    interpretation, self._not_analyte_mixture_term))
            else:
                attribs_of = interpretation.attributes

            self.handle.write(f"<Interpretation={interpretation.id}>\n")
            self._write_attributes(attribs_of)

            # When there is only one interpretation and only one interpretation member
            # interpretation member attributes are written out as part of the interpretation
            # itself.
            if _n_interps == 1 and len(interpretation.member_interpretations) == 1 and self.compact_interpretations:
                for member in interpretation.member_interpretations.values():
                    self._write_attributes(member.attributes)
            else:
                for member in interpretation.member_interpretations.values():
                    member: InterpretationMember
                    self.handle.write(f"<InterpretationMember={member.id}>\n")
                    self._write_attributes(member.attributes)
        self.handle.write("<Peaks>\n")
        for peak in spectrum.peak_list:
            peak_parts = [
                str(peak[0]),
                str(peak[1]),
                '?' if not peak[2] else ",".join(map(str, peak[2]))
            ]
            if peak[3]:
                peak_parts.append('\t'.join(map(format_aggregation, peak[3])))
            self.handle.write("\t".join(peak_parts) + "\n")
        self.handle.write("\n")

    def write_cluster(self, cluster: SpectrumCluster):
        self.handle.write(f"<Cluster={cluster.key}>\n")
        attribs_of = list(self._filter_attributes(
            cluster,
            self._not_entry_key_or_index)
        )
        self._write_attributes(attribs_of)
        self.handle.write("\n")

    def close(self):
        self.handle.close()


def format_aggregation(value: Union[float, str]) -> str:
    if isinstance(value, float):
        return "%0.4g" % value
    else:
        return value


def format_spectrum(spectrum: Spectrum, **kwargs) -> str:
    buffer = io.StringIO()
    writer = TextSpectralLibraryWriter(buffer, **kwargs)
    writer.write_spectrum(spectrum)
    return buffer.getvalue()

