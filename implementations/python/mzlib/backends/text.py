import re
import os
import io
import logging

from mzlib.index import MemoryIndex
from mzlib.annotation import parse_annotation
from mzlib.attributes import AttributeManager

from .base import (
    _PlainTextSpectralLibraryBackendBase,
    SpectralLibraryWriterBase,
    FORMAT_VERSION_TERM)
from .utils import try_cast

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


term_pattern = re.compile(
    r"^(?P<term>(?P<term_accession>\S+:\d+)\|(?P<term_name>[^=]+))")
key_value_term_pattern = re.compile(
    r"^(?P<term>(?P<term_accession>[A-Za-z0-9:.]+:\d+)\|(?P<term_name>[^=]+))=(?P<value>.+)")
grouped_key_value_term_pattern = re.compile(
    r"^\[(?P<group_id>\d+)\](?P<term>(?P<term_accession>\S+:\d+)\|(?P<term_name>[^=]+))=(?P<value>.+)")
float_number = re.compile(
    r"^\d+(.\d+)?")


START_OF_SPECTRUM_MARKER = re.compile(r"^<Spectrum>")
START_OF_ANALYTE_MARKER = re.compile(r"^<Analyte(?:=(.+))>")
START_OF_PEAKS_MARKER = re.compile(r"^<Peaks>")
START_OF_LIBRARY_MARKER = re.compile(r"^<mzSpecLib\s+(.+)>")
SPECTRUM_NAME_PRESENT = re.compile(r'MS:1003061\|spectrum name=')


class TextSpectralLibrary(_PlainTextSpectralLibraryBackendBase):
    file_format = "mzlb.txt"
    format_name = "text"

    @classmethod
    def guess_from_header(cls, filename):
        with open(filename, 'r') as stream:
            first_line = stream.readline()
            if START_OF_SPECTRUM_MARKER.match(first_line) or START_OF_LIBRARY_MARKER.match(first_line):
                return True
        return False

    def _parse_header_from_stream(self, stream):
        nbytes = 0
        first_line = stream.readline()
        nbytes += len(first_line)
        if SPECTRUM_NAME_PRESENT.match(first_line) or START_OF_SPECTRUM_MARKER.match(first_line):
            return True, 0
        elif START_OF_LIBRARY_MARKER.match(first_line):
            match = START_OF_LIBRARY_MARKER.match(first_line)
            version = match.group(1)
            attributes = AttributeManager()
            attributes.add_attribute(FORMAT_VERSION_TERM, version)
            line = stream.readline()
            while not (SPECTRUM_NAME_PRESENT.match(line) or START_OF_SPECTRUM_MARKER.match(line)):
                nbytes += len(line)
                match = key_value_term_pattern.match(line)
                if match is not None:
                    d = match.groupdict()
                    attributes.add_attribute(
                        d['term'], try_cast(d['value']))
                    line = stream.readline()
                    nbytes += len(line)
                    continue
                if line.startswith("["):
                    match = grouped_key_value_term_pattern.match(line)
                    if match is not None:
                        d = match.groupdict()
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
                    name, value = line.split("=")
                    attributes.add_attribute(name, value)
                else:
                    raise ValueError(f"Malformed attribute line {line}")
                line = stream.readline()
            self.attributes.clear()
            self.attributes._from_iterable(attributes)
            return True, nbytes
        return False, 0

    def read_header(self):
        with open(self.filename, 'rt') as stream:
            return self._parse_header_from_stream(stream)

    def create_index(self):
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

        with open(filename, 'r') as infile:
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

            logger.info(f"Reading {filename} ({file_size} bytes)...")
            while 1:
                line = infile.readline()
                if len(line) == 0:
                    break

                line_beginning_file_offset = file_offset

                #### tell() is twice as slow as counting it myself
                file_offset += len(line) + file_offset_line_ending

                line = line.rstrip()
                if state == 'header':
                    # if re.match(r'MS:1003061\|spectrum name=', line):
                    if START_OF_SPECTRUM_MARKER.match(line):
                        state = 'body'
                        spectrum_file_offset = line_beginning_file_offset
                    else:
                        continue
                if state == 'body':
                    if len(line) == 0:
                        continue
                    # if re.match(r'MS:1003061\|spectrum name=', line):
                    if START_OF_SPECTRUM_MARKER.match(line):
                        if len(spectrum_buffer) > 0:
                            if not spectrum_name:
                                raise ValueError("No spectrum name")
                            self.index.add(
                                number=n_spectra + start_index,
                                offset=spectrum_file_offset,
                                name=spectrum_name,
                                analyte=None)
                            n_spectra += 1
                            spectrum_buffer = []
                            #### Commit every now and then
                            if n_spectra % 1000 == 0:
                                self.index.commit()
                                percent_done = int(
                                    file_offset/file_size*100+0.5)
                                logger.info(str(percent_done)+"%...")

                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = ''
                    if re.match(r'MS:1003061\|spectrum name', line):
                        spectrum_name = re.match(r'MS:1003061\|spectrum name=(.+)', line).group(1)

                    spectrum_buffer.append(line)


            if not spectrum_name:
                raise ValueError("No spectrum name")
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

    def _buffer_from_stream(self, infile):
        state = 'body'
        spectrum_buffer = []

        for line in infile:
            line = line.rstrip()
            if state == 'body':
                if len(line) == 0:
                    continue
                if START_OF_SPECTRUM_MARKER.match(line):
                    if len(spectrum_buffer) > 0:
                        return spectrum_buffer
                spectrum_buffer.append(line)
        return spectrum_buffer

    def _parse(self, buffer, spectrum_index=None):
        spec = self._new_spectrum()
        analyte = None
        state = 'header'

        peak_list = []

        for line in buffer:
            line = line.strip()
            if not line:
                break
            if state == 'header':
                if START_OF_SPECTRUM_MARKER.match(line):
                    continue
                elif START_OF_PEAKS_MARKER.match(line):
                    state = 'peaks'
                    if analyte is not None:
                        spec.analytes.append(analyte)
                        analyte = None
                    continue
                elif START_OF_ANALYTE_MARKER.match(line):
                    state = 'analyte'
                    match = START_OF_ANALYTE_MARKER.match(line)
                    if analyte is not None:
                        spec.analytes.append(analyte)
                    analyte = self._new_analyte(match.group(1))
                    continue

                match = key_value_term_pattern.match(line)
                if match is not None:
                    d = match.groupdict()
                    spec.add_attribute(d['term'], try_cast(d['value']))
                    continue
                if line.startswith("["):
                    match = grouped_key_value_term_pattern.match(line)
                    if match is not None:
                        d = match.groupdict()
                        spec.add_attribute(
                            d['term'], try_cast(d['value']), d['group_id'])
                        spec.group_counter = int(d['group_id'])
                        continue
                    else:
                        raise ValueError(f"Malformed grouped attribute {line}")
                elif "=" in line:
                    name, value = line.split("=")
                    spec.add_attribute(name, value)
                else:
                    raise ValueError(f"Malformed attribute line {line}")
            elif state == 'analyte':
                if START_OF_PEAKS_MARKER.match(line):
                    state = 'peaks'
                    if analyte is not None:
                        spec.analytes.append(analyte)
                        analyte = None
                    continue
                elif START_OF_ANALYTE_MARKER.match(line):
                    state = 'analyte'
                    match = START_OF_ANALYTE_MARKER.match(line)
                    if analyte is None:
                        spec.analytes.append(analyte)
                    analyte = self._new_analyte(match.group(1))

                match = key_value_term_pattern.match(line)
                if match is not None:
                    d = match.groupdict()
                    analyte.add_attribute(d['term'], try_cast(d['value']))
                    continue
                if line.startswith("["):
                    match = grouped_key_value_term_pattern.match(line)
                    if match is not None:
                        d = match.groupdict()
                        analyte.add_attribute(
                            d['term'], try_cast(d['value']), d['group_id'])
                        analyte.group_counter = int(d['group_id'])
                        continue
                    else:
                        raise ValueError(f"Malformed grouped attribute {line}")
                elif "=" in line:
                    name, value = line.split("=")
                    analyte.add_attribute(name, value)
                else:
                    raise ValueError(f"Malformed attribute line {line}")
            elif state == 'peaks':
                match = float_number.match(line)
                if match is not None:
                    tokens = line.split("\t")
                    n_tokens = len(tokens)
                    if n_tokens == 3:
                        mz, intensity, interpretation = tokens
                        interpretation = parse_annotation(interpretation)
                        peak_list.append([float(mz), float(intensity), interpretation, ""])
                    elif n_tokens == 4:
                        mz, intensity, interpretation, aggregation = tokens
                        interpretation = parse_annotation(interpretation)
                        peak_list.append(
                            [float(mz), float(intensity), interpretation, aggregation])
                    else:
                        raise ValueError(f"Malformed peak line {line} with {n_tokens} entries")
                else:
                    raise ValueError(f"Malformed peak line {line}")
            else:
                raise ValueError(f"Unknown state {state}")
        spec.peak_list = peak_list
        return spec

    def get_spectrum(self, spectrum_number=None, spectrum_name=None):
        # keep the two branches separate for the possibility that this is not possible with all
        # index schemes.
        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError(
                    "Provide only one of spectrum_number or spectrum_name")
            offset = self.index.offset_for(spectrum_number)
        elif spectrum_name is not None:
            offset = self.index.offset_for(spectrum_name)
        buffer = self._get_lines_for(offset)
        spectrum = self._parse(buffer, spectrum_number)
        return spectrum


class TextSpectralLibraryWriter(SpectralLibraryWriterBase):
    file_format = "mzlb.txt"
    format_name = "text"
    default_version = '1.0'

    def __init__(self, filename, version=None, **kwargs):
        super(TextSpectralLibraryWriter, self).__init__(filename)
        self.version = version
        self._coerce_handle(self.filename)

    def write_header(self, library):
        if self.version is None:
            version = library.attributes.get_by_name("format version")
            if version is None:
                version = self.default_version
        else:
            version = self.version
        self.handle.write("<mzSpecLib %s>\n" % (version, ))
        for attribute in library.attributes:
            if len(attribute) == 2:
                self.handle.write(f"{attribute[0]}={attribute[1]}\n")
            elif len(attribute) == 3:
                self.handle.write(
                    f"[{attribute[2]}]{attribute[0]}={attribute[1]}\n")
            else:
                raise ValueError(
                    f"Attribute has wrong number of elements: {attribute}")

    def write_spectrum(self, spectrum):
        self.handle.write("<Spectrum>\n")
        for attribute in spectrum.attributes:
            if len(attribute) == 2:
                self.handle.write(f"{attribute[0]}={attribute[1]}\n")
            elif len(attribute) == 3:
                self.handle.write(f"[{attribute[2]}]{attribute[0]}={attribute[1]}\n")
            else:
                raise ValueError(f"Attribute has wrong number of elements: {attribute}")
        for analyte in spectrum.analytes:
            self.handle.write("<Analyte=%s>\n" % analyte.id)
            for attribute in analyte.attributes:
                if len(attribute) == 2:
                    self.handle.write(f"{attribute[0]}={attribute[1]}\n")
                elif len(attribute) == 3:
                    self.handle.write(
                        f"[{attribute[2]}]{attribute[0]}={attribute[1]}\n")
                else:
                    raise ValueError(
                        f"Attribute has wrong number of elements: {attribute}")
        self.handle.write("<Peaks>\n")
        for peak in spectrum.peak_list:
            peak_parts = [
                str(peak[0]),
                str(peak[1]),
                '?' if not peak[2] else ",".join(map(str, peak[2]))
            ]
            if peak[3]:
                peak_parts.append(str(peak[3]))
            self.handle.write("\t".join(peak_parts)+"\n")
        self.handle.write("\n")

    def close(self):
        self.handle.close()


def format_spectrum(spectrum, **kwargs):
    buffer = io.StringIO()
    writer = TextSpectralLibraryWriter(buffer, **kwargs)
    writer.write_spectrum(spectrum)
    return buffer.getvalue()

