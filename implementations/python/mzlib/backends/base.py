import os

from mzlib.index import MemoryIndex
from mzlib.spectrum import Spectrum
from mzlib.attributes import AttributeManager


class SpectralLibraryBackendBase(object):
    """A base class for all spectral library formats.

    """

    def __init__(self, filename):
        self.filename = filename
        self.index = MemoryIndex()
        self.attributes = AttributeManager()

    def add_attribute(self, key, value, group_identifier=None):
        return self.attributes.add_attribute(key, value, group_identifier=group_identifier)

    def get_attribute(self, key, group_identifier=None):
        return self.attributes.get_attribute(key, group_identifier=group_identifier)

    def _new_spectrum(self):
        return Spectrum()

    def get_spectrum(self, spectrum_number=None, spectrum_name=None):
        raise NotImplementedError()

    def search(self, specification, **query_keys):
        raise NotImplementedError()

    def create_index(self):
        """
        Populate the spectrum index.

        This method may produce a large amount of file I/O.

        Returns
        -------
        n_spectra: int
            The number of entries read
        """
        raise NotImplementedError()

    def __iter__(self):
        for record in self.index:
            yield self.get_spectrum(record.number)

    def __len__(self):
        return len(self.index)

    def __getitem__(self, i):
        record = self.index[i]
        if isinstance(record, list):
            result = [self.get_spectrum(rec.number) for rec in record]
        else:
            result = self.get_spectrum(record.number)
        return result


class _PlainTextSpectralLibraryBackendBase(SpectralLibraryBackendBase):

    def __init__(self, filename, index_type=None, read_metadata=True):
        if index_type is None:
            index_type = MemoryIndex
        super(_PlainTextSpectralLibraryBackendBase, self).__init__(filename)
        self.index, was_initialized = index_type.from_filename(filename)
        if not was_initialized:
            self.create_index()
        if read_metadata:
            self.read_header()

    def read_header(self):
        raise NotImplementedError()

    def _get_lines_for(self, offset):
        raise NotImplementedError()

    def _parse(self, buffer, spectrum_index=None):
        raise NotImplementedError()

    def search(self, specification, **query_keys):
        records = self.index.search(specification, **query_keys)
        if not isinstance(records, list):
            records = [records]
        spectra = []
        for record in records:
            buffer = self._get_lines_for(record.offset)
            spectrum = self._parse(buffer, record.number)
            spectra.append(spectrum)
        return spectra
