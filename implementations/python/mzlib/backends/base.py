import os

from mzlib.index import MemoryIndex
from mzlib.spectrum import Spectrum


class SpectralLibraryBackendBase(object):

    def __init__(self, filename):
        self.filename = filename

    def _new_spectrum(self):
        return Spectrum()

    def get_spectrum(self, spectrum_number=None, spectrum_name=None):
        raise NotImplementedError()

    def search(self, specification, **query_keys):
        raise NotImplementedError()


class _PlainTextSpectralLibraryBackendBase(SpectralLibraryBackendBase):

    def __init__(self, filename, index_type=None):
        if index_type is None:
            index_type = MemoryIndex
        super(_PlainTextSpectralLibraryBackendBase, self).__init__(filename)
        self.index, was_initialized = index_type.from_filename(filename)
        if not was_initialized:
            self._build_index()

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
