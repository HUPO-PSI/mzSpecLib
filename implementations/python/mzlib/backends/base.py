import os

from mzlib.index import MemoryIndex
from mzlib.spectrum import Spectrum
from mzlib.analyte import Analyte
from mzlib.attributes import AttributeManager


class SpectralLibraryBackendBase(object):
    """A base class for all spectral library formats.

    """

    def __init__(self, filename):
        self.filename = filename
        self.index = MemoryIndex()
        self.attributes = AttributeManager()

    def read_header(self):
        """read_header - Read just the header of the whole library

        Returns
        -------
        bool
        """
        raise NotImplementedError()

    def add_attribute(self, key, value, group_identifier=None):
        """Add an attribute to the library level attributes store.

        Parameters
        ----------
        key : str
            The name of the attribute to add
        value : object
            The value of the attribute to add
        group_identifier : str, optional
            The attribute group identifier to use, if any. If not provided,
            no group is assumed.
        """
        return self.attributes.add_attribute(key, value, group_identifier=group_identifier)

    def get_attribute(self, key, group_identifier=None):
        """Get the value or values associated with a given
        attribute key from the library level attribute store.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        Returns
        -------
        attribute_value: object or list[object]
            Returns single or multiple values for the requested attribute.
        """
        return self.attributes.get_attribute(key, group_identifier=group_identifier)

    def remove_attribute(self, key, group_identifier=None):
        """Remove the value or values associated with a given
        attribute key from the library level attribute store.

        This rebuilds the entire store, which may be expensive.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        """
        return self.attributes.remove_attribute(key, group_identifier=group_identifier)

    def has_attribute(self, key):
        """Test for the presence of a given attribute in the library
        level store.

        Parameters
        ----------
        key : str
            The attribute to test for

        Returns
        -------
        bool
        """
        return self.attributes.has_attribute(key)

    def _new_spectrum(self):
        return Spectrum()

    def _new_analyte(self, id=None):
        return Analyte(id)

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
        raise NotImplementedError()

    def find_spectra(self, specification, **query_keys):
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

    def _coerce_handle(self, filename_or_stream):
        if hasattr(filename_or_stream, 'read'):
            self.handle = filename_or_stream
        else:
            self.handle = open(filename_or_stream, 'rt')

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


class SpectralLibraryWriterBase(object):
    def __init__(self, filename, **kwargs):
        self.filename = filename

    def _coerce_handle(self, filename_or_stream):
        if hasattr(filename_or_stream, 'write'):
            self.handle = filename_or_stream
        else:
            self.handle = open(filename_or_stream, 'wt')

    def write_library(self, library):
        self.write_header(library)
        for spectrum in library:
            self.write_spectrum(spectrum)

    def write_spectrum(self, spectrum):
        raise NotImplementedError()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        pass
