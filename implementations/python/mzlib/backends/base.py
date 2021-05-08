import os

from pathlib import Path

from mzlib.index import MemoryIndex, SQLIndex
from mzlib.spectrum import Spectrum
from mzlib.analyte import Analyte
from mzlib.attributes import AttributeManager


FORMAT_VERSION_TERM = 'MS:1009002|format version'
DEFAULT_VERSION = '1.0'


class SubclassRegisteringMetaclass(type):
    def __new__(mcs, name, parents, attrs):
        new_type = type.__new__(mcs, name, parents, attrs)
        if not hasattr(new_type, "_file_extension_to_implementation"):
            new_type._file_extension_to_implementation = dict()

        file_extension = attrs.get("file_format")
        if file_extension is not None:
            new_type._file_extension_to_implementation[file_extension] = new_type

        format_name = attrs.get("format_name")
        if format_name is not None:
            new_type._file_extension_to_implementation[format_name] = new_type
        else:
            attrs['format_name'] = file_extension
        return new_type

    def type_for_format(cls, format_or_extension):
        return cls._file_extension_to_implementation.get(format_or_extension)


class SpectralLibraryBackendBase(object, metaclass=SubclassRegisteringMetaclass):
    """A base class for all spectral library formats.

    """
    file_format = None

    _file_extension_to_implementation = {}

    @classmethod
    def guess_from_filename(cls, filename):
        """Guess if the file is of this type by inspecting the file's name and extension.

        Parameters
        ----------
        filename : str
            The path to the file to inspect.

        Returns
        -------
        bool:
            Whether this is an appropriate backend for that file.
        """
        if not isinstance(filename, (str, Path)):
            return False
        return filename.endswith(cls.file_format)

    @classmethod
    def guess_from_header(cls, filename):
        """Guess if the file is of this type by inspecting the file's header section

        Parameters
        ----------
        filename : str
            The path to the file to open.

        Returns
        -------
        bool:
            Whether this is an appropriate backend for that file.
        """
        return False

    @classmethod
    def guess_implementation(cls, filename, index_type=None, **kwargs):
        """Guess the backend implementation to use with this file format.

        Parameters
        ----------
        filename : str
            The path to the spectral library file to open.
        index_type : type, optional
            The :class:`~.IndexBase` derived type to use for this file. If
            :const:`None` is provided, the instance will decide based upon
            :meth:`has_index_preference`.

        Returns
        -------
        SpectralLibraryBackendBase
        """
        for key, impl in cls._file_extension_to_implementation.items():
            try:
                if impl.guess_from_filename(filename):
                    return impl(filename, index_type=index_type, **kwargs)
            except TypeError:
                pass
            try:
                if impl.guess_from_header(filename):
                    return impl(filename, index_type=index_type, **kwargs)
            except TypeError:
                pass
        raise ValueError(f"Could not guess backend implementation for {filename}")

    def __init__(self, filename):
        self.filename = filename
        self.index = MemoryIndex()
        self.attributes = AttributeManager()

    @property
    def format_version(self):
        try:
            value = self.get_attribute(FORMAT_VERSION_TERM)
            return value
        except KeyError:
            value = DEFAULT_VERSION
            self.add_attribute(FORMAT_VERSION_TERM, value)
            return value

    def read_header(self):
        """Read just the header of the whole library

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
        """Populate the spectrum index.

        This method may produce a large amount of file I/O.

        Returns
        -------
        n_spectra: int
            The number of entries read
        """
        raise NotImplementedError()

    def __iter__(self):
        if self.index:
            for record in self.index:
                yield self.get_spectrum(record.number)
        else:
            raise NotImplementedError()
            return self.read()

    def __len__(self):
        return len(self.index)

    def __getitem__(self, i):
        record = self.index[i]
        if isinstance(record, list):
            result = [self.get_spectrum(rec.number) for rec in record]
        else:
            result = self.get_spectrum(record.number)
        return result

    @classmethod
    def has_index_preference(cls, filename):
        '''Does this backend prefer a particular index for this file?

        The base implementation checks to see if there is a SQL index
        for the filename provided, and if so, prefers :class:`~.SQLIndex`.
        Otherwise, prefers :class:`~.MemoryIndex`.

        Parameters
        ----------
        filename: str
            The name of the file to open.

        Returns
        -------
        index_type: type
            Returns a :class:`~.IndexBase` derived type which this backend
            would prefer to use.
        '''
        try:
            if SQLIndex.exists(filename):
                return SQLIndex
            return MemoryIndex
        except Exception:
            return MemoryIndex

    def read(self):
        raise NotImplementedError()

guess_implementation = SpectralLibraryBackendBase.guess_implementation


class _PlainTextSpectralLibraryBackendBase(SpectralLibraryBackendBase):

    def __init__(self, filename, index_type=None, read_metadata=True):
        if index_type is None:
            index_type = self.has_index_preference(filename)
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

    def _buffer_from_stream(self, stream):
        '''Collect data from the readable stream until
        a complete spectrum entry has been observed.

        Parameters
        ----------
        stream: file-like
            Theinput file stream to read from.

        Returns
        -------
        line_buffer: list[str]
            A list of lines read from the input stream.
        '''
        raise NotImplementedError()

    def read(self):
        with open(self.filename, 'rt') as stream:
            i = 0
            match, offset = self._parse_header_from_stream(stream)
            if not match:
                raise ValueError("Could not locate valid header")
            else:
                stream.seek(offset)
            while True:
                # Will clip the first line of the next spectrum. Needs work
                buffer = self._buffer_from_stream(stream)
                if not buffer:
                    break
                yield self._parse(buffer, i)

    def _get_lines_for(self, offset):
        with open(self.filename, 'r') as infile:
            infile.seek(offset)
            spectrum_buffer = self._buffer_from_stream(infile)
            #### We will end up here if this is the last spectrum in the file
        return spectrum_buffer

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


class SpectralLibraryWriterBase(object, metaclass=SubclassRegisteringMetaclass):
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
