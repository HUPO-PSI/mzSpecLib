import io
import enum
import logging

from typing import Callable, Dict, Iterable, Union, List, Type, Iterator
from pathlib import Path


from psims.controlled_vocabulary import Entity
from psims.controlled_vocabulary.controlled_vocabulary import (
    load_uo, load_unimod, load_psims)

from mzlib.index import MemoryIndex, SQLIndex, IndexBase
from mzlib.spectrum import LIBRARY_ENTRY_INDEX, LIBRARY_ENTRY_KEY, Spectrum
from mzlib.analyte import Analyte, Interpretation, InterpretationMember, ANALYTE_MIXTURE_TERM
from mzlib.attributes import Attributed, AttributedEntity, AttributeSet, AttributeManagedProperty
from mzlib.ontology import _VocabularyResolverMixin

from .utils import open_stream, LineBuffer

logger = logging.getLogger(__name__.rsplit(".", 1)[0])
logger.addHandler(logging.NullHandler())

ANALYTE_MIXTURE_CURIE = ANALYTE_MIXTURE_TERM.split("|")[0]

FORMAT_VERSION_TERM = 'MS:1003186|library format version'
LIBRARY_NAME_TERM = "MS:1003188|library name"
LIBRARY_VERSION_TERM = "MS:1003190|library version"
LIBRARY_IDENTIFIER_TERM = "MS:1003187|library identifier"
LIBRARY_DESCRIPTION_TERM = "MS:1003189|library description"
LIBRARY_URI_TERM = "MS:1003191|library URI"


DEFAULT_VERSION = '1.0'


class AttributeSetTypes(enum.Enum):
    spectrum = enum.auto()
    analyte = enum.auto()
    interpretation = enum.auto()
    cluster = enum.auto()


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


class SpectralLibraryBackendBase(AttributedEntity, _VocabularyResolverMixin, metaclass=SubclassRegisteringMetaclass):
    """A base class for all spectral library formats.

    """
    file_format = None

    _file_extension_to_implementation: Dict[str,
                                            Type['SpectralLibraryBackendBase']] = {}
    _format_name_to_implementation: Dict[str,
                                         Type['SpectralLibraryBackendBase']] = {}

    index: IndexBase

    entry_attribute_sets: Dict[str, AttributeSet]
    analyte_attribute_sets: Dict[str, AttributeSet]
    interpretation_attribute_sets: Dict[str, AttributeSet]
    cluster_attribute_sets: Dict[str, AttributeSet]

    name = AttributeManagedProperty[str](LIBRARY_NAME_TERM)
    identifier = AttributeManagedProperty[str](LIBRARY_IDENTIFIER_TERM)
    description = AttributeManagedProperty[str](LIBRARY_DESCRIPTION_TERM)
    uri = AttributeManagedProperty[str](LIBRARY_URI_TERM)

    @classmethod
    def guess_from_filename(cls, filename: Union[str, Path, io.FileIO]) -> bool:
        """
        Guess if the file is of this type by inspecting the file's name and extension.

        Parameters
        ----------
        filename : str
            The path to the file to inspect.

        Returns
        -------
        bool:
            Whether this is an appropriate backend for that file.
        """
        if hasattr(filename, "name"):
            filename = filename.name
        if not isinstance(filename, (str, Path)):
            return False
        if filename.endswith(".gz"):
            filename = filename[:-3]
        return filename.endswith(cls.file_format)

    @classmethod
    def guess_from_header(cls, filename) -> bool:
        """
        Guess if the file is of this type by inspecting the file's header section

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
    def guess_implementation(cls, filename, index_type=None,
                             **kwargs) -> 'SpectralLibraryBackendBase':
        """
        Guess the backend implementation to use with this file format.

        Parameters
        ----------
        filename : str
            The path to the spectral library file to open.
        index_type : type, optional
            The :class:`~.IndexBase` derived type to use for this file. If
            :const:`None` is provided, the instance will decide based upon
            :meth:`has_index_preference`.
        **kwargs
            Passed to implementation

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
            except (TypeError, UnicodeDecodeError):
                pass
        raise ValueError(f"Could not guess backend implementation for {filename}")

    def __init__(self, filename):
        self.filename = filename
        self.index = MemoryIndex()

        self.entry_attribute_sets = {
            "all": AttributeSet("all", [])
        }
        self.analyte_attribute_sets = {
            "all": AttributeSet("all", [])
        }
        self.interpretation_attribute_sets = {
            "all": AttributeSet("all", [])
        }

        super().__init__(None)

    @property
    def format_version(self):
        try:
            value = self.get_attribute(FORMAT_VERSION_TERM)
            return value
        except KeyError:
            value = DEFAULT_VERSION
            self.add_attribute(FORMAT_VERSION_TERM, value)
            return value

    def read_header(self) -> bool:
        """
        Read just the header of the whole library

        Returns
        -------
        bool
        """
        raise NotImplementedError()

    def _new_spectrum(self) -> Spectrum:
        spec = Spectrum()
        attr_set = self.entry_attribute_sets.get("all")
        if attr_set:
            attr_set.apply(spec)
        return spec

    def _new_interpretation(self, id=None) -> Interpretation:
        interp = Interpretation(id)
        attr_set = self.interpretation_attribute_sets.get('all')
        if attr_set:
            attr_set.apply(interp)
        return interp

    def _new_interpretation_member(self, id=None) -> InterpretationMember:
        return InterpretationMember(id)

    def _new_analyte(self, id=None) -> Analyte:
        analyte = Analyte(id)
        attr_set = self.analyte_attribute_sets.get('all')
        if attr_set:
            attr_set.apply(analyte)
        return analyte

    def _new_cluster(self) -> SpectrumCluster:
        cluster = SpectrumCluster()
        attr_set = self.cluster_attribute_sets.get('all')
        if attr_set:
            attr_set.apply(cluster)
        return cluster

    def _analyte_interpretation_link(self, spectrum: Spectrum,
                                     interpretation: Interpretation):
        if (interpretation.has_attribute(ANALYTE_MIXTURE_TERM) and
            not interpretation.analytes):
            analyte_ids = interpretation.get_attribute(ANALYTE_MIXTURE_TERM)
            if isinstance(analyte_ids, str):
                term = self.find_term_for(ANALYTE_MIXTURE_CURIE)
                analyte_ids = term.value_type(analyte_ids)

            # TODO: Enforce this attribute is a string at the CV level
            # if isinstance(analyte_ids_term, int):
            #     analyte_ids = [analyte_ids_term]
            #     interpretation.replace_attribute(ANALYTE_MIXTURE_TERM, str(analyte_ids_term))
            # else:
            #     analyte_ids = analyte_ids_term.split(',')
            for analyte_id in analyte_ids:
                interpretation.add_analyte(spectrum.get_analyte(analyte_id))
        return interpretation

    def _default_interpretation_to_analytes(self, spectrum: Spectrum):
        for interpretation in spectrum.interpretations.values():
            if not interpretation.analytes:
                for analyte in spectrum.analytes.values():
                    interpretation.add_analyte(analyte)

    def get_spectrum(self, spectrum_number: int=None,
                     spectrum_name: str=None) -> Spectrum:
        """
        Retrieve a single spectrum from the library.

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

    def get_cluster(self, cluster_number: int) -> SpectrumCluster:
        raise NotImplementedError()

    def find_spectra(self, specification, **query_keys):
        raise NotImplementedError()

    def create_index(self) -> int:
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
        if self.index:
            for record in self.index:
                yield self.get_spectrum(record.number)
        else:
            return self.read()

    def __len__(self):
        return len(self.index)

    def __getitem__(self, i) -> Union[Spectrum, List[Spectrum]]:
        record = self.index[i]
        if isinstance(record, list):
            result = [self.get_spectrum(rec.number) for rec in record]
        else:
            result = self.get_spectrum(record.number)
        return result

    @classmethod
    def has_index_preference(cls, filename: str) -> Type[IndexBase]:
        """
        Does this backend prefer a particular index for this file?

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
        """
        try:
            if SQLIndex.exists(filename):
                return SQLIndex
            return MemoryIndex
        except Exception:
            return MemoryIndex

    def read(self):
        raise NotImplementedError()

    def _add_attribute_set(self, attribute_set: AttributeSet,
                           attribute_set_type: AttributeSetTypes):
        if attribute_set_type == AttributeSetTypes.spectrum:
            self.entry_attribute_sets[attribute_set.name] = attribute_set
        elif attribute_set_type == AttributeSetTypes.analyte:
            self.analyte_attribute_sets[attribute_set.name] = attribute_set
        elif attribute_set_type == AttributeSetTypes.interpretation:
            self.interpretation_attribute_sets[attribute_set.name] = attribute_set
        elif attribute_set_type == AttributeSetTypes.cluster:
            self.cluster_attribute_sets[attribute_set.name] = attribute_set
        else:
            raise ValueError(f"Could not map {attribute_set_type}")

    def summarize_parsing_errors(self) -> Dict:
        return {}

guess_implementation = SpectralLibraryBackendBase.guess_implementation


class _PlainTextSpectralLibraryBackendBase(SpectralLibraryBackendBase):

    def __init__(self, filename, index_type=None, read_metadata=True,
                 create_index: bool=True):
        if index_type is None and create_index:
            index_type = self.has_index_preference(filename)

        super(_PlainTextSpectralLibraryBackendBase, self).__init__(filename)

        if index_type is not None:
            self.index, was_initialized = index_type.from_filename(filename)
            if not was_initialized and create_index:
                self.create_index()
        if read_metadata:
            self.read_header()

    def _coerce_handle(self, filename_or_stream):
        if hasattr(filename_or_stream, 'read'):
            self.handle = filename_or_stream
        else:
            self.handle = open_stream(filename_or_stream, 'rt')

    def _buffer_from_stream(self, stream: io.IOBase) -> List:
        """
        Collect data from the readable stream until
        a complete spectrum entry has been observed.

        Parameters
        ----------
        stream: file-like
            Theinput file stream to read from.

        Returns
        -------
        line_buffer: List[str]
            A list of lines read from the input stream.
        """
        raise NotImplementedError()

    def read(self) -> Iterator[Spectrum]:
        with open_stream(self.filename, 'rt') as stream:
            i = 0
            match, offset = self._parse_header_from_stream(stream)
            if not match:
                raise ValueError("Could not locate valid header")
            else:
                stream.seek(offset)
            buffering_stream = LineBuffer(stream)
            while True:
                # Will clip the first line of the next spectrum. Needs work
                buffer = self._buffer_from_stream(buffering_stream)

                # If the buffer is only a single line, then we must have reached
                # the end, so we're done. We're done because the buffering stream
                # will contain exactly one line (the buffered line)
                if len(buffer) <= 1:
                    break
                buffering_stream.push_line()
                yield self._parse(buffer, i)
                i += 1

    def _get_lines_for(self, offset: int) -> List[str]:
        filename = self.filename
        is_file_like_object = isinstance(filename, io.IOBase)

        infile = open_stream(filename, 'r')

        infile.seek(offset)
        spectrum_buffer = self._buffer_from_stream(infile)
        #### We will end up here if this is the last spectrum in the file
        if not is_file_like_object:
            infile.close()
        else:
            infile.detach()
        return spectrum_buffer

    def _parse(self, buffer: Iterable, spectrum_index: int=None):
        raise NotImplementedError()

    def search(self, specification, **query_keys) -> List[Spectrum]:
        records = self.index.search(specification, **query_keys)
        if not isinstance(records, list):
            records = [records]
        spectra = []
        for record in records:
            buffer = self._get_lines_for(record.offset)
            spectrum = self._parse(buffer, record.number)
            spectra.append(spectrum)
        return spectra


class SpectralLibraryWriterBase(_VocabularyResolverMixin, metaclass=SubclassRegisteringMetaclass):
    def __init__(self, filename, **kwargs):
        self.filename = filename
        super().__init__(**kwargs)

    def _filter_attributes(self, attributes: Attributed,
                           filter_fn: Callable) -> Iterable:
        if isinstance(attributes, AttributedEntity):
            attributes = attributes.attributes
        for attrib in attributes:
            if filter_fn(attrib):
                yield attrib

    def _not_analyte_mixture_term(self, attrib):
        if attrib:
            key = attrib[0]
            if key == ANALYTE_MIXTURE_TERM:
                return False
        return True

    def _not_entry_index(self, attrib):
        if attrib:
            key = attrib[0]
            if key == LIBRARY_ENTRY_INDEX:
                return False
        return True

    def _not_entry_key_or_index(self, attrib):
        if attrib:
            key = attrib[0]
            if key in (LIBRARY_ENTRY_INDEX, LIBRARY_ENTRY_KEY):
                return False
        return True

    def _coerce_handle(self, filename_or_stream):
        if hasattr(filename_or_stream, 'write'):
            self.handle = filename_or_stream
        else:
            self.handle = open(filename_or_stream, 'wt')

    def write_library(self, library: SpectralLibraryBackendBase):
        self.write_header(library)
        n = len(library)
        step = max(min(n // 100, 5000), 1)
        ident = ''
        i = 0
        for i, spectrum in enumerate(library):
            if i % step == 0 and i:
                try:
                    ident = f"{spectrum.key}:{spectrum.name}"
                except Exception:
                    ident = str(spectrum.key)
                logger.info(f"Wrote {ident} {i}/{n} ({i / n * 100.0:0.2f}%)")
            self.write_spectrum(spectrum)
        i = n
        logger.info(f"Wrote {n} spectra")

    def write_spectrum(self, spectrum: Spectrum):
        raise NotImplementedError()

    def write_cluster(self, cluster: SpectrumCluster):
        raise NotImplementedError()

    def __enter__(self) -> 'SpectralLibraryWriterBase':
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        pass


class LibrarySpectrumIterator(AttributedEntity, Iterator[Spectrum]):
    def __init__(self, backend: SpectralLibraryBackendBase) -> None:
        self.backend = backend
        self.attributes = backend
        self.iter = backend.read()
        self._buffer = next(self.iter)

    @property
    def format_version(self):
        return self.backend.format_version

    def __iter__(self):
        return self

    def __next__(self) -> Spectrum:
        if self._buffer is not None:
            result = self._buffer
            self._buffer = None
            return result
        return next(self.iter)
