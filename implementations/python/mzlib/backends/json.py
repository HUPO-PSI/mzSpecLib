import io
import json
import re

from pathlib import Path

from mzlib.index import MemoryIndex
from mzlib.attributes import AttributeManager
from mzlib.annotation import parse_annotation, IonAnnotationBase

from .base import SpectralLibraryBackendBase, SpectralLibraryWriterBase, FORMAT_VERSION_TERM


LIBRARY_METADATA_KEY = "attributes"
ELEMENT_ATTRIBUTES_KEY = "attributes"
LIBRARY_SPECTRA_KEY = "spectra"
FORMAT_VERSION_KEY = "format_version"

FORMAT_VERSION_ACC = FORMAT_VERSION_TERM.split("|")[0]


class JSONSpectralLibrary(SpectralLibraryBackendBase):
    file_format = "mzlb.json"
    format_name = "json"

    def __init__(self, filename, index_type=None, read_metadata=True):
        if index_type is None:
            index_type = MemoryIndex
        super(JSONSpectralLibrary, self).__init__(filename)
        self.buffer = {}
        self._load_buffer(self.filename)
        self.attributes = AttributeManager(
            self.buffer.get(LIBRARY_METADATA_KEY))
        self.index, was_initialized = index_type.from_filename(self.filename)
        if not was_initialized:
            self.create_index()

    @classmethod
    def guess_from_filename(cls, filename):
        if isinstance(filename, dict):
            return LIBRARY_SPECTRA_KEY in filename and LIBRARY_METADATA_KEY in filename
        if not isinstance(filename, (str, Path)):
            return False
        return filename.endswith(cls.file_format)

    def _load_buffer(self, filename_or_stream):
        if isinstance(filename_or_stream, dict):
            self.buffer = filename_or_stream
        else:
            if hasattr(filename_or_stream, 'read'):
                self.handle = filename_or_stream
            else:
                self.handle = open(filename_or_stream, 'rt')
            self.buffer = json.load(self.handle)
            self.handle.close()

    def read_header(self):
        if self.buffer:
            pass
        return False

    def create_index(self):
        for i, record in enumerate(self.buffer[LIBRARY_SPECTRA_KEY]):
            for attrib in record['attributes']:
                if attrib["accession"] == "MS:1003061":
                    self.index.add(i, i, attrib['value'], None, None)
                    break
            else:
                raise ValueError(f"Unidentified spectrum at index {i}")

    def get_spectrum(self, spectrum_number=None, spectrum_name=None):
        """Retrieve a single spectrum from the library.

        Parameters
        ----------
        spectrum_number : int, optional
            The index of the specturm in the library
        spectrum_name : str, optional
            The name of the spectrum in the library

        Returns
        -------
        :class:`~.Spectrum`
        """

        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError(
                    "Provide only one of spectrum_number or spectrum_name")
            offset = self.index.offset_for(spectrum_number)
        elif spectrum_name is not None:
            offset = self.index.offset_for(spectrum_name)
        data = self.buffer[LIBRARY_SPECTRA_KEY][offset]
        spectrum = self.make_spectrum_from_payload(data)
        return spectrum

    def make_spectrum_from_payload(self, data):
        spectrum = self._new_spectrum()
        for attrib in data['attributes']:
            key = f'{attrib["accession"]}|{attrib["name"]}'
            if "value_accession" in attrib:
                value = f'{attrib["value_accession"]}|{attrib["value"]}'
            else:
                value = attrib['value']
            group = attrib.get("cv_param_group")
            spectrum.add_attribute(key, value, group_identifier=group)
            if group is not None:
                spectrum.group_counter = int(group)
        analytes = []
        for analyte_id, analyte in data['analytes'].items():
            analyte_d = self._new_analyte(analyte_id)
            for attrib in analyte['attributes']:
                key = f'{attrib["accession"]}|{attrib["name"]}'
                if "value_accession" in attrib:
                    value = f'{attrib["value_accession"]}|{attrib["value"]}'
                else:
                    value = attrib['value']
                group = attrib.get("cv_param_group")
                analyte_d.add_attribute(key, value, group_identifier=group)
                if group is not None:
                    analyte_d.group_counter = int(group)
            analytes.append(analyte_d)
        spectrum.analytes = analytes
        peak_list = []
        n = len(data['mzs'])
        mzs = data['mzs']
        intensities = data['intensities']
        interpretations = data['interpretations']
        aggregations = data.get("aggregations", None)
        for i in range(n):
            interpretation = interpretations[i]
            if isinstance(interpretation, str):
                interpretation = parse_annotation(interpretation)
            elif isinstance(interpretation, list):
                interpretation = [IonAnnotationBase.from_json(interp) for interp in interpretation]
            elif isinstance(interpretation, dict):
                interpretation = [IonAnnotationBase.from_json(interpretation)]
            else:
                raise TypeError(f"Cannot reconstruct interpretation from type {interpretation.__class__}")
            peak = [
                mzs[i],
                intensities[i],
                parse_annotation(interpretations[i]),
                aggregations[i] if aggregations else ''
            ]
            peak_list.append(peak)
        spectrum.peak_list = peak_list
        return spectrum

    def read(self):
        n = len(self.buffer[LIBRARY_SPECTRA_KEY])
        for offset in range(n):
            data = self.buffer[LIBRARY_SPECTRA_KEY][offset]
            spectrum = self.make_spectrum_from_payload(data)
            yield spectrum


class JSONSpectralLibraryWriter(SpectralLibraryWriterBase):
    file_format = "mzlb.json"
    format_name = "json"
    default_version = '1.0'

    def __init__(self, filename, version=None, pretty_print=True, format_annotations=True, simplify=True, **kwargs):
        if version is None:
            version = self.default_version
        super(JSONSpectralLibraryWriter, self).__init__(filename)
        self._coerce_handle(self.filename)
        self.version = version
        self.pretty_print = pretty_print
        self.wrote_library = False
        self.simplify = simplify
        self.format_annotations = format_annotations
        self.buffer = {
            FORMAT_VERSION_KEY: self.version,
            LIBRARY_METADATA_KEY: [],
            LIBRARY_SPECTRA_KEY: []
        }

    def write_library(self, library):
        self.wrote_library = True
        return super().write_library(library)

    def split_compound_value(self, value):
        value = str(value)
        # Don't process quoted values
        if value.startswith('"'):
            return [value]
        components = value.split('|', 1)
        return components

    def write_header(self, library):
        attributes = []
        for attribute in library.attributes:
            reformed_attribute = {}
            if len(attribute) == 2:
                key, value = attribute
            elif len(attribute) == 3:
                key, value, cv_param_group = attribute
                reformed_attribute['cv_param_group'] = cv_param_group
            else:
                raise ValueError(
                    f"Unsupported number of items in attribute: {attribute}")
            if FORMAT_VERSION_ACC == key:
                # Already covered by format_version
                continue
            components = key.split('|', 1)
            if len(components) == 2:
                accession,name = components
                reformed_attribute['accession'] = accession
                reformed_attribute['name'] = name
            else:
                raise ValueError(
                    f"Unsupported number of items in components: {components}")

            components = self.split_compound_value(value)
            if len(components) == 2:
                value_accession,value = components
                reformed_attribute['value_accession'] = value_accession
                reformed_attribute['value'] = value
            elif len(components) == 1:
                reformed_attribute['value'] = value
            else:
                raise ValueError(
                    f"Unsupported number of items in components: {components}")
            attributes.append(reformed_attribute)
        self.buffer[LIBRARY_METADATA_KEY] = attributes

    def _format_attributes(self, attributes_manager):
        attributes = []
        for attribute in attributes_manager:
            reformed_attribute = {}
            if len(attribute) == 2:
                key, value = attribute
            elif len(attribute) == 3:
                key, value, cv_param_group = attribute
                reformed_attribute['cv_param_group'] = cv_param_group
            else:
                raise ValueError(
                    f"Unsupported number of items in attribute: {attribute}")

            components = key.split('|', 1)
            if len(components) == 2:
                accession, name = components
                reformed_attribute['accession'] = accession
                reformed_attribute['name'] = name
            else:
                raise ValueError(
                    f"Unsupported number of items in components: {components}")

            components = self.split_compound_value(value)
            if len(components) == 2:
                value_accession, value = components
                reformed_attribute['value_accession'] = value_accession
                reformed_attribute['value'] = value
            elif len(components) == 1:
                reformed_attribute['value'] = value
            else:
                raise ValueError(
                    f"Unsupported number of items in components: {components}")
            attributes.append(reformed_attribute)
        return attributes

    def write_spectrum(self, spectrum):
        mzs = []
        intensities = []
        interpretations = []
        aggregations = []
        for peak in spectrum.peak_list:
            mzs.append(peak[0])
            intensities.append(peak[1])
            if self.format_annotations:
                interpretations.append(
                    '?' if not peak[2] else ",".join(map(str, peak[2])))
            else:
                interpretations.append([c.to_json() for c in peak[2]])
            aggregations.append(peak[3])

        #### Organize the attributes from the simple list into the appropriate JSON format
        attributes = self._format_attributes(spectrum)

        analytes = {}
        for analyte in spectrum.analytes:
            analyte_d = {
                "id": analyte.id,
                ELEMENT_ATTRIBUTES_KEY: self._format_attributes(analyte)
            }
            analytes[analyte.id] = (analyte_d)

        spectrum = {
            ELEMENT_ATTRIBUTES_KEY: attributes,
            "mzs": mzs,
            "intensities": intensities,
            "interpretations": interpretations,
            "aggregations": aggregations,
            "analytes": analytes
        }
        if not any(aggregations):
            spectrum.pop('aggregations')

        self.buffer[LIBRARY_SPECTRA_KEY].append(spectrum)

    def flush(self):
        # If we know we're writing a complete library, skip the probably-doing-too-many-things
        # formatting logic for single vs. many spectra.
        if self.wrote_library:
            if self.pretty_print:
                json.dump(self.buffer, self.handle, indent=2, sort_keys=True)
            else:
                json.dump(self.buffer, self.handle)
        else:
            # We don't have a header section to format, so write just the spectra,
            # and if the number of spectra is one and the simplify flag is true,
            # skip the wrapping array
            spectra = self.buffer[LIBRARY_SPECTRA_KEY]
            n_spectra = len(spectra)
            if n_spectra == 1 and self.simplify:
                if self.pretty_print:
                    json.dump(spectra[0], self.handle,
                              indent=2, sort_keys=True)
                else:
                    json.dump(spectra[0], self.handle)
            else:
                if self.pretty_print:
                    json.dump(spectra, self.handle, indent=2, sort_keys=True)
                else:
                    json.dump(spectra, self.handle)

    def close(self):
        self.flush()
        self.handle.close()


def format_spectrum(spectrum, pretty_print=True, **kwargs):
    buffer = io.StringIO()
    with JSONSpectralLibraryWriter(buffer, pretty_print=pretty_print, **kwargs) as writer:
        writer.write_spectrum(spectrum)
        writer.flush()
        return buffer.getvalue()
